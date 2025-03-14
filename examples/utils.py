import streamlit as st
import numpy as np
from simcc import *

materialOptions = [
    "Be",
    "Al",
    "Mylar",
    "Pla",
    "Kapton",
    "P10",
    "P10 620 Torr",
    "Xe7",
    "Xe7 620 Torr",
    "CH4",
    "CH4 620 Torr",
    "Diamond",
    "Gold",
    "Tungsten",
]

# 各カテゴリの物質リスト（厚みの範囲を考慮）
material_list = {
    "Gas detectors": ["P10 gas IC", "Xe gas IC", "PPAC"],
    "Gas": ["P10", "Xe7", "iC4H10", "CH4", "GasHe", "GasNe", "GasAr", "GasKr", "GasXe"],
    "Strippers": ["Al", "Ta", "W", "Pt", "Au"],
    "Degraders": ["Al", "Cu", "Acrylic"],
    "Targets": ["Be", "W", "Carbon"],
    "Others": ["Plastic", "Diamond", "Kapton", "Mylar"],
    "Single substance": [f"Z={Z}" for Z in range(1, 93)],
}

clight = 299.792458  # mm/ns
mnucleon = 931.49432  # MeV


def energy2beta(E):
    KE = E
    M = mnucleon
    P = ((M + KE) ** 2 - M**2) ** 0.5
    gamma2 = 1 + (P / mnucleon) ** 2
    beta = (1 - (1 / gamma2)) ** 0.5
    return beta


def beta2energy(beta):
    gamma2 = 1 / (1 - beta * beta)
    P = mnucleon * (gamma2 - 1) ** 0.5
    M = mnucleon
    KE = (P * P + M * M) ** 0.5 - M
    E = KE
    return E


def beta2mom(beta, A):
    gamma = (1 / (1 - beta * beta)) ** 0.5
    return beta * gamma * A * mnucleon


def beta2brho(beta, A, Q):
    return beta2mom(beta, A) / Q / clight


def brho2beta(brho, A, Q):
    P = brho * Q * clight
    gamma2 = 1 + (P / (A * mnucleon)) ** 2
    beta = (1 - (1 / gamma2)) ** 0.5
    return beta


def brho2energy(brho, A, Q):
    return beta2energy(brho2beta(brho, A, Q))


def energy2brho(energy, A, Q):
    return beta2brho(energy2beta(energy), A, Q)


def input_projectile():
    st.write("**Projectile**: Configure the projectile parameters as the initial state for the calculations.")
    use_number_input = st.toggle("Number input")

    col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
    with col1:
        if use_number_input:
            projectile_Z = st.number_input("Z", min_value=30, max_value=94, step=1, value=70)
        else:
            projectile_Z = st.slider("Z", 30, 94, 70)
        st.write(f"Element: {z2symbol[projectile_Z]}")
    with col2:
        if use_number_input:
            energy = st.number_input("Energy (MeV/u)", min_value=50.0, max_value=1000.0, step=5.0, value=250.0)
        else:
            energy = st.slider("Energy (MeV/u)", 50, 1000, 250, step=5)
    with col3:
        if use_number_input:
            A = st.number_input("A", min_value=50, max_value=300, step=1, value=175)
        else:
            AoZ = st.slider("A/Z", 1.5, 3.5, 2.5, step=0.01)
            A = int(np.round(projectile_Z * AoZ))
            st.write(f"A={A}")
    with col4:
        charge_states = {0: "Full-strip", 1: "H-like", 2: "He-like", 3: "Li-like", 4: "Be-like", 5: "B-like", 6: "C-like"}
        charge_state = st.selectbox("Charge states", options=list(charge_states.keys()), format_func=lambda x: charge_states[x])
        st.write(f"Brho={energy2brho(energy, A, projectile_Z-charge_state):.6f} Tm")

    return projectile_Z, energy, A, charge_state


def input_materials():
    st.write("**Material Setting**: Select and configure materials. Press 'Add' to append them to the material list.")

    material_units = ["mm", "μm", "cm", "mg/cm2"]
    material_unit = st.radio("Unit", material_units, key=f"material_unit", horizontal=True)

    # 選択された物質のリスト (セッションステートを使用)
    if "selected_materials" not in st.session_state:
        st.session_state.selected_materials = ["P10 gas IC-0"]
    if "j" not in st.session_state:
        st.session_state.j = 1

    # 横に並べる
    col1, col2, col3, col4 = st.columns([3, 2, 2, 1])

    with col1:
        category = st.selectbox("Category", list(material_list.keys()))

    with col2:
        material = st.selectbox("Material", material_list[category])
        if category != "Gas detectors":
            density = GetMaterial(material)["density"]
            if category == "Gas":
                st.write(f"{density*1000:.6g} mg/cm3")
            elif category == "Single substance":
                st.write(f"Element: {z2symbol[int(material.replace("Z=", ""))]} ({density:.6g} g/cm3)")
            elif category != "Gas detectors":
                st.write(f"{density:.6g} g/cm3")

    with col3:
        if category == "Gas detectors":
            thickness = st.number_input("Thickness (N/A)", disabled=True)  # 入力不可
            expanded_materials = get_expanded_materials([material])
            st.write(", ".join(expanded_materials))
        elif material_unit == "mg/cm2":
            thickness_mgcm3 = st.number_input("Thickness (mg/cm2)", min_value=0.0, max_value=5000.0, step=1.0, value=100.0)
            st.session_state.thickness = thickness_mgcm3 / density / 1000 / 0.1
        else:
            if material_unit == "μm":
                thickness = st.number_input("Thickness (µm)", min_value=0.0, max_value=1000.0, step=10.0, value=10.0)
                st.session_state.thickness = thickness * 0.001
            elif material_unit == "mm":
                thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=500.0, step=1.0, value=1.0)
                st.session_state.thickness = thickness * 1
            elif material_unit == "cm":
                thickness = st.number_input("Thickness (cm)", min_value=0.1, max_value=500.0, step=1.0, value=1.0)
                st.session_state.thickness = thickness * 10
            if category == "Gas":
                st.write(f"{density*1000 * st.session_state.thickness*0.1:.6g} mg/cm2")
            else:
                st.write(f"{density * st.session_state.thickness*0.1:.6g} g/cm2")

    with col4:
        if st.button("Add"):
            if category == "Gas detectors":
                item = f"{material}"
            elif category == "Single substance":
                item = f"{material} {density * st.session_state.thickness * 0.1:.6g} g/cm2"
            else:
                item = f"{material} {st.session_state.thickness:.6g} mm"
            st.session_state.selected_materials.append(f"{item}-{st.session_state.j}")
            st.session_state.j += 1

    st.write("**Material List**: These are used in the calculations. Uncheck to remove.")

    for i, item in enumerate(st.session_state.selected_materials):
        checked = st.checkbox(item.split("-")[0], value=True, key=item)
        if checked == False:
            st.session_state.selected_materials.pop(i)
            st.rerun()
    return st.session_state.selected_materials


def get_expanded_materials(materials):
    expanded_materials = []
    for i, material in enumerate(materials):
        material = material.split("-")[0]
        if material in material_list["Gas detectors"]:
            if material == material_list["Gas detectors"][0]:
                expanded_materials.append("Kapton 0.125 mm")
                expanded_materials.append("P10 586 mm")
                expanded_materials.append("Mylar 0.1 mm")
                expanded_materials.append("Kapton 0.125 mm")
            elif material == material_list["Gas detectors"][1]:
                expanded_materials.append("Kapton 0.125 mm")
                expanded_materials.append("Xe7 586 mm")
                expanded_materials.append("Mylar 0.1 mm")
                expanded_materials.append("Kapton 0.125 mm")
            elif material == material_list["Gas detectors"][2]:
                expanded_materials.append("Mylar 0.045 mm")
        else:
            expanded_materials.append(material)
    return expanded_materials


def get_material_name_length(material):
    material_name = material.split()[0]
    if material.split()[2] == "mm":
        length = float(material.split()[1])
    elif material.split()[2] == "g/cm2":
        length = float(material.split()[1])*10
    return material_name, length
