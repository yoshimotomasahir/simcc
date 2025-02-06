import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
import sys

sys.path.append("..")
from simcc import *
from utils import *

st.write("Projectile:")
col2, col3 = st.columns(2)
with col2:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=70, step=1)
with col3:
    energy = st.number_input("Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=250.0, step=5.0)

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    EqDist = []
    MFPs = []
    Zs = np.arange(30, 93)
    for Zt in Zs:
        MFP = GetMFP(Z, energy, int(Zt), solid_gas="solid")
        MFPs.append(MFP)
        EqDist.append(GetEqDist(MFP))

    fig, ax = plt.subplots(figsize=(5, 2))
    for i, val in enumerate(np.array(EqDist).T):
        ax.plot(Zs, val, label=f"Q={Z-i}")
    ax.set_title("Equilibrium charge distribution")
    ax.set_xlabel("Target atomic number Z")
    ax.grid(alpha=0.4)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    st.pyplot(fig)

    for label1, label2 in zip([f"{Z}->{Z-1}", f"{Z-1}->{Z-2}", f"{Z-2}->{Z-3}"], [f"{Z-1}->{Z}", f"{Z-2}->{Z-1}", f"{Z-3}->{Z-2}"]):
        fig, ax = plt.subplots(figsize=(5, 2))
        MFP_array = np.array([np.array(list(MFP.values())) for MFP in MFPs]).T
        ax.plot(Zs, [MFP[label1] * 1000 for MFP in MFPs], label=label1)
        ax.plot(Zs, [MFP[label2] * 1000 for MFP in MFPs], label=label2)
        ax.set_title("Mean free path [mg/cm$^2$]")
        ax.set_xlabel("Target atomic number Z")
        ax.grid(alpha=0.4)
        # ax.set_ylim(0, None)
        ax.set_yscale("log")
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        st.pyplot(fig)
