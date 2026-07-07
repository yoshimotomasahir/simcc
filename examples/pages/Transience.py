import streamlit as st
import numpy as np
from simcc import *
import sys

sys.path.append("..")
from utils import *

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd

st.set_page_config(page_title="Transient  state -SimCC-", page_icon="🌠")
st.header("Transient  state")

projectile = input_projectile(use_url_params=True)
projectile_Z = projectile["Z"]
energy = projectile["energy"]
A = projectile["A"]
mass = projectile["mass"]
charge_state = projectile["charge_state"]

exp_correction = input_exp_correction()

materials = input_materials()

if st.button("Execute Calculation"):

    P0 = np.zeros(7)
    P0[charge_state] = 1
    normalized_lengths = [0]
    actual_lengths = [0]
    Probs = [P0]
    z_effectives = [projectile_Z, projectile_Z - 1, projectile_Z - 2, projectile_Z - 3, 1]
    energies = {}

    for z_effective in z_effectives:
        energies[z_effective] = [energy]

    expanded_materials = get_expanded_materials(materials)
    material_end_indices = []
    material_labels = []

    n = len(expanded_materials)
    for i, material in enumerate(expanded_materials):
        material_name, length = get_material_name_length(material)

        energy0 = {z_effective: energies[z_effective][-1] for z_effective in z_effectives}

        Eloss = GetAnalyticalEloss(A, projectile_Z, energy0[1], material_name, length * 0.1)[0]
        length_log = np.array([0] + list(np.logspace(np.log10(length / 1000), np.log10(length), num=20 + int(Eloss))))

        for l0, l1 in zip(length_log[:-1], length_log[1:]):
            normalized_lengths.append(i + l1 / length)
            actual_lengths.append(actual_lengths[-1] + l1 - l0)

            for z_effective in z_effectives:
                energy1 = energy0[z_effective] - GetAnalyticalEloss(A, projectile_Z, energy0[z_effective], material_name, l1 * 0.1, z_effective=z_effective)[0]
                if z_effective == 1:
                    MFP = GetMFP(zp=projectile_Z, energy=energy1, material=material_name, exp_correction=exp_correction)
                    P0 = GetAnalyticalProb(MFP, (l1 - l0) * 0.1, charge_state=P0)
                    Probs.append(P0)
                energies[z_effective].append(energy1)
        material_end_indices.append(len(normalized_lengths) - 1)
        material_labels.append(material.split("-")[0])

    # グラフ・テーブル化
    probs = np.array(Probs)

    fig = make_subplots(subplot_titles=["Charge-state probability"])
    for j in range(7):
        fig.add_trace(go.Scatter(x=normalized_lengths, y=probs.T[j].T, mode="lines", name=f"Q={projectile_Z-j}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=1, text=expanded_materials[i].split("-")[0], showarrow=False)
    st.plotly_chart(fig)

    table_indices = [0] + material_end_indices
    table_labels = ["Initial"] + [f"{i + 1}: {material_label}" for i, material_label in enumerate(material_labels)]

    prob_rows = []
    for j in range(7):
        row = {"Charge state": f"q={projectile_Z-j}"}
        for table_label, end_index in zip(table_labels, table_indices):
            row[table_label] = probs[end_index][j] * 100
        prob_rows.append(row)
    st.dataframe(pd.DataFrame(prob_rows).set_index("Charge state").round(3), width="stretch")

    fig = make_subplots(subplot_titles=["Energy (MeV/u) by CATIMA"])
    for z_effective in z_effectives:
        fig.add_trace(go.Scatter(x=normalized_lengths, y=energies[z_effective], mode="lines", name=f"q=Zeff" if z_effective == 1 else f"q={z_effective}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=energy, text=expanded_materials[i].split("-")[0], showarrow=False)
    st.plotly_chart(fig)

    energy_rows = []
    for z_effective in z_effectives:
        label = "q=Zeff" if z_effective == 1 else f"q={z_effective}"
        row = {"Charge state": f"{label}"}
        for table_label, end_index in zip(table_labels, table_indices):
            row[table_label] = energies[z_effective][end_index]
        energy_rows.append(row)
    st.dataframe(pd.DataFrame(energy_rows).set_index("Charge state").round(4), width="stretch")

    header = ["Actual_length", "Normalized_length"] + ["Energy_Zeff"] + [f"q={projectile_Z-j}" for j in range(7)]
    data_rows = [
        [f"{actual_lengths[j]:.6g}"] + [f"{normalized_lengths[j]:.6g}"] + [f"{energies[1][j]:.5f}"] + [f"{p:.6e}" for p in probs[j]]
        for j in range(len(normalized_lengths))
    ]
    st.text_area("Raw data", "\n".join(["\t".join(data) for data in [header] + data_rows]), height=300)

    st.success("Calculation executed successfully!")
