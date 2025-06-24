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
    "Others": ["Plastic", "Diamond", "Kapton", "Mylar", "EmulsionE07"],
    "Single substance": [f"Z={Z}" for Z in range(1, 93)],
}

F5_deg_list = [
    "Nothing 0 mm 0 mrad",
    "Al 0.7 mm 0.61 mrad",
    "Al 1.0 mm 0.43 mrad",
    "Al 1.0 mm 0.95 mrad",
    "Al 1.5 mm 1.187 mrad",
    "Al 2.0 mm 1.60 mrad",
    "Al 2.2 mm 1.80 mrad",
    "Al 2.5 mm 2.14 mrad",
    "Al 3.0 mm 2.62 mrad",
    "Al 3.5 mm 2.80 mrad",
    "Al 3.5 mm 3.00 mrad",
    "Al 4.5 mm 3.70 mrad",
    "Al 5.0 mm 4.25 mrad",
    "Al 7.0 mm 5.969 mrad",
    "Al 8.0 mm 7.3129 mrad",
    "Al 9.0 mm 7.714 mrad",
    "Al 10.0 mm 8.587 mrad",
]

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


def tof2beta(FL, TOF):
    beta = FL / TOF / clight
    return beta


def beta2tof(FL, beta):
    TOF = FL / clight / beta
    return TOF


def get_matrix(Fa, Fb):
    import requests

    url = "https://ribf.riken.jp/BigRIPSInfo/optics/fig/bStandard.txt"
    response = requests.get(url)
    assert response.status_code == 200
    data = response.text
    lines = data.splitlines()
    mat = []
    for i, line in enumerate(lines):
        if line.replace("\n", "").replace("\r", "") == "F{}-F{}".format(Fa, Fb):
            for j in range(6):
                mat.append([float(f) for f in lines[i + j + 1].split()[0:5]])
            assert lines[i + 7][0] == "-"
            return mat


def input_projectile(comment=""):
    st.write("**Projectile**: Configure the projectile parameters as the initial state for the calculations.")
    if comment != "":
        st.write(comment)
    use_number_input = st.toggle("Number input", value=True)

    col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
    with col2:
        if use_number_input:
            projectile_Z = st.number_input("Z", min_value=30, max_value=94, step=1, value=70)
        else:
            projectile_Z = st.slider("Z", 30, 94, 70)
        st.write(f"Element: {z2symbol[projectile_Z]}")
    with col4:
        energy_unit = st.radio("Energy unit", ["MeV/u", "Tm"], horizontal=True)
        if use_number_input:
            if energy_unit == "MeV/u":
                energy = st.number_input("Energy (MeV/u)", min_value=50.0, max_value=1000.0, step=5.0, value=300.0)
            elif energy_unit == "Tm":
                brho = st.number_input("Energy (Tm)", min_value=2.0, max_value=20.0, step=0.1, value=6.7)
        else:
            if energy_unit == "MeV/u":
                energy = st.slider("Energy (MeV/u)", 50, 1000, 300, step=5)
            elif energy_unit == "Tm":
                brho = st.number_input("Energy (Tm)", min_value=2.0, max_value=20.0, step=0.1, value=6.7)
    with col1:
        if use_number_input:
            A = st.number_input("A", min_value=50, max_value=300, step=1, value=175)
        else:
            AoZ = st.slider("A/Z", 1.5, 3.5, 2.5, step=0.01)
            A = int(np.round(projectile_Z * AoZ))
            st.write(f"A={A}")
    with col3:
        charge_states = {0: "Full-strip", 1: "H-like", 2: "He-like", 3: "Li-like", 4: "Be-like", 5: "B-like", 6: "C-like"}
        charge_state = st.selectbox("Charge states", options=list(charge_states.keys()), format_func=lambda x: charge_states[x])
        if energy_unit == "MeV/u":
            st.write(f"{energy2brho(energy, A, projectile_Z-charge_state):.6f} Tm")
        else:
            energy = brho2energy(brho, A, projectile_Z - charge_state)
            st.write(f"{energy:.6f} MeV/u")

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


def input_exp_correction():
    cs_option = st.sidebar.selectbox("Version of charge state changing cross section", ["Original v0", "Correction v1"])
    exp_correction = int(cs_option.split()[1].replace("v", ""))
    return exp_correction


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
        length = float(material.split()[1]) * 10
    return material_name, length
