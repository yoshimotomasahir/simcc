import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
from simcc import *
import sys

sys.path.append("..")
from utils import *

import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

st.title("Equilibrium Thickness")

st.write("The charge-state equilibrium is considered achieved when the difference from the equilibrium state probability is less than 1/e^6.")

graph_type = st.selectbox("Select the graph to display", ["Projectile Z Dependence", "Energy Dependence", "Target Z Dependence"])

solid_gas = "solid"

def get_user_inputs(graph_type):

    col1, col2, col3 = st.columns(3)

    with col1:
        if graph_type == "Projectile Z Dependence":
            projectile_Z = st.slider("Projectile Z Range", 30, 92, (30, 92))
        else:
            projectile_Z = st.slider("Projectile Z", 30, 92, 70)
    with col2:
        if graph_type == "Energy Dependence":
            energy = st.slider("Energy (MeV/u) Range", 50, 1000, (100, 350))
        else:
            energy = st.slider("Energy (MeV/u)", 50, 1000, 250)
    with col3:
        if graph_type == "Target Z Dependence":
            target_Z = st.slider("Target Z Range", 1, 92, (1, 92))
        else:
            target_Z = st.slider("Target Z", 1, 92, 13)
    return projectile_Z, energy, target_Z

projectile_Z, energy, target_Z = get_user_inputs(graph_type)

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    if graph_type == "Energy Dependence":
        energies = np.arange(energy[0], energy[1], 10)
        thickness = {}
        for charge_state in range(3):
            thickness[charge_state] = [GetEquilibriumThickness(GetMFP(projectile_Z, ene, target_Z, solid_gas=solid_gas), charge_state) for ene in energies]
        EqDist = [GetEqDist(GetMFP(projectile_Z, ene, target_Z, solid_gas="solid")) for ene in energies]

        fig = make_subplots(subplot_titles=["Equilibrium Thickness (mg/cm2)"])
        for charge_state in range(3):
            fig.add_trace(go.Scatter(x=energies, y=np.array(thickness[charge_state]) * 1000, name=f"Qin = {projectile_Z-charge_state}"))
        fig.update_layout(yaxis=dict(range=(0, max(np.array([t for t in thickness.values()]).flatten()) * 1000 * 1.1)))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)

        fig = make_subplots(subplot_titles=["Charge-state distribution"])
        for i, val in enumerate(np.array(EqDist).T):
            fig.add_trace(go.Scatter(x=energies, y=val, mode='lines', name=f"Q={projectile_Z-i}"))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)

    elif graph_type == "Projectile Z Dependence":
        pzs = np.arange(projectile_Z[0], projectile_Z[1] + 1)
        thickness = {}
        for charge_state in range(3):
            thickness[charge_state] = [GetEquilibriumThickness(GetMFP(pz, energy, target_Z, solid_gas=solid_gas), charge_state) for pz in pzs]
        EqDist = [GetEqDist(GetMFP(pz, energy, target_Z, solid_gas="solid")) for pz in pzs]

        fig = make_subplots(subplot_titles=["Equilibrium Thickness (mg/cm2)"])
        for charge_state in range(3):
            fig.add_trace(go.Scatter(x=pzs, y=np.array(thickness[charge_state]) * 1000, name=f"Qin-Z = {charge_state}"))
        fig.update_layout(yaxis=dict(range=(0, max(np.array([t for t in thickness.values()]).flatten()) * 1000 * 1.1)))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)

        fig = make_subplots(subplot_titles=["Charge-state distribution"])
        for i, val in enumerate(np.array(EqDist).T):
            fig.add_trace(go.Scatter(x=pzs, y=val, mode='lines', name=f"Z-Q={i}"))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)

    elif graph_type == "Target Z Dependence":
        tzs = np.arange(target_Z[0], target_Z[1] + 1).tolist()
        thickness = {}
        for charge_state in range(3):
            thickness[charge_state] = [GetEquilibriumThickness(GetMFP(projectile_Z, energy, tz, solid_gas=solid_gas), charge_state) for tz in tzs]
        EqDist = [GetEqDist(GetMFP(projectile_Z, energy, tz, solid_gas="solid")) for tz in tzs]

        fig = make_subplots(subplot_titles=["Equilibrium Thickness (mg/cm2)"])
        for charge_state in range(3):
            fig.add_trace(go.Scatter(x=tzs, y=np.array(thickness[charge_state]) * 1000, name=f"Qin = {projectile_Z-charge_state}"))
        fig.update_layout(yaxis=dict(range=(0, max(np.array([t for t in thickness.values()]).flatten()) * 1000 * 1.1)))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)

        fig = make_subplots(subplot_titles=["Charge-state distribution"])
        for i, val in enumerate(np.array(EqDist).T):
            fig.add_trace(go.Scatter(x=tzs, y=val, mode='lines', name=f"Z-Q={i}"))
        fig.update_layout(xaxis_title=" ".join(graph_type.split()[:-1]))
        st.plotly_chart(fig)
