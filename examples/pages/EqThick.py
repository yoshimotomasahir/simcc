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

st.title("Equilibrium Thickness")

st.write("The charge-state equilibrium is considered achieved when the difference from the equilibrium state probability is less than 1/e^6.")

graph_type = st.selectbox("Select the graph to display", ["Projectile Z Dependence", "Energy Dependence", "Target Z Dependence"])

charge_state = 0
solid_gas = "solid"

def get_user_inputs(graph_type):

    if graph_type == "Projectile Z Dependence":
        projectile_Z = st.slider("Projectile Z Range", 30, 92, (30, 92))
    else:
        projectile_Z = st.slider("Projectile Z", 30, 92, 70)

    if graph_type == "Energy Dependence":
        energy = st.slider("Energy (MeV/u) Range", 50, 1000, (100, 350))
    else:
        energy = st.slider("Energy (MeV/u)", 50, 1000, 250)

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
        thickness = [GetEquilibriumThickness(GetMFP(projectile_Z, ene, target_Z, solid_gas=solid_gas), charge_state) for ene in energies]
        fig = px.line(x=energies, y=np.array(thickness) * 1000, labels={"x": "Energy (MeV)", "y": "Equilibrium Thickness (mg/cm2)"}, title="Energy Dependence")
        fig.update_layout(yaxis=dict(range=(0, max(thickness) * 1000 * 1.1)))

    elif graph_type == "Projectile Z Dependence":
        pzs = np.arange(projectile_Z[0], projectile_Z[1] + 1)
        thickness = [GetEquilibriumThickness(GetMFP(pz, energy, target_Z, solid_gas=solid_gas), charge_state) for pz in pzs]
        fig = px.line(x=pzs, y=np.array(thickness) * 1000, labels={"x": "Projectile Z", "y": "Equilibrium Thickness (mg/cm2)"}, title="Projectile Z Dependence")

    elif graph_type == "Target Z Dependence":
        tzs = np.arange(target_Z[0], target_Z[1] + 1).tolist()
        thickness = [GetEquilibriumThickness(GetMFP(projectile_Z, energy, tz, solid_gas=solid_gas), charge_state) for tz in tzs]
        fig = px.line(x=tzs, y=np.array(thickness) * 1000, labels={"x": "Target Z", "y": "Equilibrium Thickness (mg/cm2)"}, title="Target Z Dependence")

    st.plotly_chart(fig)
