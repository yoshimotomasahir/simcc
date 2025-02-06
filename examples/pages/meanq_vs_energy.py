import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
import sys

sys.path.append("..")
from simcc import *
from utils import *

st.write("Projectile:")
col1, col2, col3 , col4 = st.columns(4)
with col1:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=70, step=1)
with col2:
    energy0 = st.number_input("Min Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=150.0, step=5.0)
with col3:
    energy1 = st.number_input("Max Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=350.0, step=5.0)
with col4:
    energyBin = st.number_input("Number of bins", min_value=2, max_value=100, value=10, step=1)

st.write("Target:")
col1, = st.columns(1)
with col1:
    Zt = st.number_input("Target atomic number Z", min_value=1, max_value=120, value=73, step=1)

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    EqDist = []
    MFPs = []
    energies = np.linspace(energy0, energy1, energyBin)
    for energy in energies:
        MFP = GetMFP(Z, energy, int(Zt), solid_gas="solid")
        MFPs.append(MFP)
        EqDist.append(GetEqDist(MFP))

    fig, ax = plt.subplots(figsize=(5, 2))
    for i, val in enumerate(np.array(EqDist).T):
        ax.plot(energies, val, label=f"Q={Z-i}")
    ax.set_title("Equilibrium charge distribution")
    ax.set_xlabel("Energy [MeV/u]")
    ax.grid(alpha=0.4)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    st.pyplot(fig)

    for label1, label2 in zip([f"{Z}->{Z-1}", f"{Z-1}->{Z-2}", f"{Z-2}->{Z-3}"], [f"{Z-1}->{Z}", f"{Z-2}->{Z-1}", f"{Z-3}->{Z-2}"]):
        fig, ax = plt.subplots(figsize=(5, 2))
        MFP_array = np.array([np.array(list(MFP.values())) for MFP in MFPs]).T
        ax.plot(energies, [MFP[label1] * 1000 for MFP in MFPs], label=label1)
        ax.plot(energies, [MFP[label2] * 1000 for MFP in MFPs], label=label2)
        ax.set_title("Mean free path [mg/cm$^2$]")
        ax.set_xlabel("Energy [MeV/u]")
        ax.grid(alpha=0.4)
        # ax.set_ylim(0, None)
        ax.set_yscale("log")
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        st.pyplot(fig)
