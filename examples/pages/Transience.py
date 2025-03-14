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

st.set_page_config(page_title="Transient  state -SimCC-", page_icon="🌠")
st.title("Transient  state")

projectile_Z, energy, A, charge_state = input_projectile()

materials = input_materials()

if st.button("Execute Calculation"):

    P0 = np.zeros(7)
    P0[charge_state] = 1
    lengths = [0]
    Probs = [P0]
    z_effectives = [projectile_Z, projectile_Z - 1, projectile_Z - 2, projectile_Z - 3, 1]
    energies = {}

    for z_effective in z_effectives:
        energies[z_effective] = [energy]

    expanded_materials = get_expanded_materials(materials)

    n = len(expanded_materials)
    for i, material in enumerate(expanded_materials):
        material_name, length = get_material_name_length(material)

        energy0 = {z_effective: energies[z_effective][-1] for z_effective in z_effectives}

        Eloss = GetAnalyticalEloss(A, projectile_Z, energy0[1], material_name, length * 0.1)[0]
        length_log = np.array([0] + list(np.logspace(np.log10(length / 1000), np.log10(length), num=20 + int(Eloss))))

        sub_energies = []
        for l0, l1 in zip(length_log[:-1], length_log[1:]):
            lengths.append(i + l1 / length)

            for z_effective in z_effectives:
                energy1 = energy0[z_effective] - GetAnalyticalEloss(A, projectile_Z, energy0[z_effective], material_name, l1 * 0.1, z_effective=z_effective)[0]
                if z_effective == 1:
                    MFP = GetMFP(zp=projectile_Z, energy=energy1, material=material_name)
                    P0 = GetAnalyticalProb(MFP, (l1 - l0) * 0.1, charge_state=P0)
                    Probs.append(P0)
                energies[z_effective].append(energy1)

    # グラフ化
    fig = make_subplots(subplot_titles=["Energy (MeV/u) by CATIMA"])
    for z_effective in z_effectives:
        fig.add_trace(go.Scatter(x=lengths, y=energies[z_effective], mode="lines", name=f"Q=Zeff" if z_effective == 1 else f"Q={z_effective}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=energy, text=expanded_materials[i].split("-")[0], showarrow=False)
    st.plotly_chart(fig)

    fig = make_subplots(subplot_titles=["Charge-state probability"])
    for j in range(7):
        fig.add_trace(go.Scatter(x=lengths, y=np.array(Probs).T[j].T, mode="lines", name=f"Z-Q={j}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=1, text=expanded_materials[i].split("-")[0], showarrow=False)
    st.plotly_chart(fig)
    header = ["Length"] + ["Energy"] + [f"Z-Q={j}" for j in range(7)]
    data_rows = [[str(lengths[j])] + [str(energies[1][j])] + list(map(str, np.array(Probs)[j].T)) for j in range(len(lengths))]
    st.text_area("Raw data", "\n".join(["\t".join(data) for data in [header] + data_rows]), height=300)

    st.success("Calculation executed successfully!")
