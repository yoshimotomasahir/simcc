import sys
sys.path.append("..")
from utils import *

from simcc import *
import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.title("Charge-state distribution")

st.write("Projectile:")
col1, col2, col3 , col4 = st.columns(4)
with col1:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=70, step=1)
with col2:
    energy0 = st.number_input("Min Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=150.0, step=5.0)
with col3:
    energy1 = st.number_input("Max Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=350.0, step=5.0)
with col4:
    energyBin = st.number_input("Number of bins", min_value=2, max_value=100, value=100, step=1)

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

    fig = go.Figure()

    for i, val in enumerate(np.array(EqDist).T):
        fig.add_trace(go.Scatter(x=energies, y=val, mode='lines', name=f"Q={Z-i}"))

    fig.update_layout(
        title="Equilibrium charge distribution",
        xaxis_title="Energy [MeV/u]",
        yaxis_title="Equilibrium Distribution",
        template="plotly",
        legend_title="Charge States"
    )

    st.plotly_chart(fig)
