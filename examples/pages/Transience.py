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

st.write("**Projectile**")

col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
with col1:
    projectile_Z = st.slider("Projectile Z", 30, 92, 70)
with col2:
    energy = st.slider("Energy [MeV/u]", 50, 1000, 250, step=5)
with col3:
    Aoq = st.slider("A/q", 1.5, 3.5, 2.5, step=0.02)
with col4:
    charge_states = {0: "Full-strip", 1: "H-like", 2: "He-like"}
    charge_state = st.selectbox("Charge states", options=list(charge_states.keys()), format_func=lambda x: charge_states[x])

st.write("**Material**")

# 各カテゴリの物質リスト（厚みの範囲を考慮）
material_list = {
    "Strippers": ["Ta", "Al", "W", "Pt", "Au"],
    "Gas detectors": ["P10 gas IC", "Xe gas IC", "PPAC"],
    "Other detectors": ["Plastic", "Diamond"],
    "Degraders": ["Al", "Cu"],
    "Targets": ["Be", "W"],
}

# 選択された物質のリスト (セッションステートを使用)
if "selected_materials" not in st.session_state:
    st.session_state.selected_materials = []
if "j" not in st.session_state:
    st.session_state.j = 0


# 横に並べる
col1, col2, col3, col4 = st.columns([3, 2, 2, 1])

with col1:
    category = st.selectbox("Category", list(material_list.keys()))

with col2:
    material = st.selectbox("Material", material_list[category])

with col3:
    if category == "Gas detectors":
        thickness = st.number_input("Thickness (N/A)", disabled=True)  # 入力不可
    elif category == "Strippers":
        thickness = st.number_input("Thickness (µm)", min_value=0.0, max_value=1000.0, step=10.0, value=10.0)
    elif category == "Targets":
        thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=20.0, step=1.0, value=2.0)
    elif category == "Degraders":
        thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=20.0, step=0.1, value=2.0)
    elif category == "Other detectors":
        thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=10.0, step=0.1, value=0.1)

with col4:
    if st.button("Add"):
        if category == "Gas detectors":
            item = f"{material}"
        else:
            item = f"{material} {thickness*0.001 if category =='Strippers' else thickness} mm"
        st.session_state.selected_materials.append(f"{item}-{st.session_state.j}")
        st.session_state.j += 1

for i, item in enumerate(st.session_state.selected_materials):
    checked = st.checkbox(item.split("-")[0], value=True, key=item)
    if checked == False:
        st.session_state.selected_materials.pop(i)
        st.rerun()

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

    materials = []
    items = []
    for i, item in enumerate(st.session_state.selected_materials):
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
        materials.append(item.split("-")[0])
        length = float(item.split()[1])

        energy0 = {z_effective: energies[z_effective][-1] for z_effective in z_effectives}

        Eloss = GetAnalyticalEloss(int(projectile_Z * 2.5), projectile_Z, energy0[1], material, length * 0.1)[0]
        length_log = np.array([0] + list(np.logspace(np.log10(length / 1000), np.log10(length), num=20 + int(Eloss))))

        sub_energies = []
        for l0, l1 in zip(length_log[:-1], length_log[1:]):
            lengths.append(i + l1 / length)

            for z_effective in z_effectives:
                energy1 = energy0[z_effective] - GetAnalyticalEloss(int(projectile_Z * 2.5), projectile_Z, energy0[z_effective], material, l1 * 0.1, z_effective=z_effective)[0]
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
        fig.add_annotation(x=i + 0.5, y=energy, text=materials[i], showarrow=False)
    st.plotly_chart(fig)

    fig = make_subplots(subplot_titles=["Charge-state probability"])
    for j in range(4):
        fig.add_trace(go.Scatter(x=lengths, y=np.array(Probs).T[j].T, mode="lines", name=f"Z-Q={j}"))
    for j in range(n + 1):
        fig.add_vline(x=j, line_dash="dash", line_color="gray", line_width=1)
    for i in range(n):
        fig.add_annotation(x=i + 0.5, y=1, text=materials[i], showarrow=False)
    st.plotly_chart(fig)

    st.success("Calculation executed successfully!")
