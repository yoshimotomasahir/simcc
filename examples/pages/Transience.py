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

st.title("Transient  state")

projectile_Z, energy, A, charge_state = input_projectile()

materials = input_materials()

if st.button("Execute Calculation"):

    P0 = np.zeros(7)
    P0[charge_state] = 1
    energy0 = energy
    lengths = [0]
    Probs = [P0]
    z_effectives = [projectile_Z, projectile_Z - 1, projectile_Z - 2, projectile_Z - 3, 1]
    energies = {}

    for z_effective in z_effectives:
        energies[z_effective] = [energy]

    raw_materials = []
    items = []
    for i, item in enumerate(materials):
        item = item.split("-")[0]
        if item in material_list["Gas detectors"]:
            if item == material_list["Gas detectors"][0]:
                items.append("Kapton 0.125 mm")
                items.append("P10 586 mm")
                items.append("Mylar 0.1 mm")
                items.append("Kapton 0.125 mm")
            elif item == material_list["Gas detectors"][1]:
                items.append("Kapton 0.125 mm")
                items.append("Xe7 586 mm")
                items.append("Mylar 0.1 mm")
                items.append("Kapton 0.125 mm")
            elif item == material_list["Gas detectors"][2]:
                items.append("Mylar 0.045 mm")
        else:
            items.append(item)
    n = len(items)
    for i, item in enumerate(items):
        material = item.split()[0]
        raw_materials.append(item.split("-")[0])
        length = float(item.split()[1])

        energy0 = {z_effective: energies[z_effective][-1] for z_effective in z_effectives}

        Eloss = GetAnalyticalEloss(A, projectile_Z, energy0[1], material, length * 0.1)[0]
        length_log = np.array([0] + list(np.logspace(np.log10(length / 1000), np.log10(length), num=20 + int(Eloss))))

        sub_energies = []
        for l0, l1 in zip(length_log[:-1], length_log[1:]):
            lengths.append(i + l1 / length)

            for z_effective in z_effectives:
                energy1 = energy0[z_effective] - GetAnalyticalEloss(A, projectile_Z, energy0[z_effective], material, l1 * 0.1, z_effective=z_effective)[0]
                if z_effective == 1:
                    MFP = GetMFP(zp=projectile_Z, energy=energy1, material=material)
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
        fig.add_annotation(x=i + 0.5, y=energy, text=raw_materials[i], showarrow=False)
    st.plotly_chart(fig)

    fig = make_subplots(subplot_titles=["Charge-state probability"])
    for j in range(7):
        fig.add_trace(go.Scatter(x=lengths, y=np.array(Probs).T[j].T, mode="lines", name=f"Z-Q={j}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=1, text=raw_materials[i], showarrow=False)
    st.plotly_chart(fig)

    st.success("Calculation executed successfully!")
