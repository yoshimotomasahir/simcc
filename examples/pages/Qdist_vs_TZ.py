import sys
sys.path.append("..")
from utils import *

from simcc import *
import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.title("Charge-state distribution")

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

    fig = go.Figure()
    
    for i, val in enumerate(np.array(EqDist).T):
        fig.add_trace(go.Scatter(x=Zs, y=val, mode='lines', name=f"Q={Z-i}"))
    
    fig.update_layout(
        title="Equilibrium charge distribution",
        xaxis_title="Target Z",
        yaxis_title="Equilibrium Distribution",
        template="plotly",
        legend_title="Charge States"
    )
    
    st.plotly_chart(fig)