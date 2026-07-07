import streamlit as st
import numpy as np
from simcc import *
from functools import lru_cache

mnucleon = 931.49410242  # // 統一原子質量単位 MeV/c^2
emass = 0.51099895000  #   // 電子質量 MeV/c^2

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


def beta2mom(beta, mass):
    gamma = (1 / (1 - beta * beta)) ** 0.5
    return beta * gamma * mass * mnucleon


def beta2brho(beta, mass, Q):
    return beta2mom(beta, mass) / Q / clight


def brho2beta(brho, mass, Q):
    P = brho * Q * clight
    gamma2 = 1 + (P / (mass * mnucleon)) ** 2
    beta = (1 - (1 / gamma2)) ** 0.5
    return beta


def brho2energy(brho, mass, Q):
    return beta2energy(brho2beta(brho, mass, Q))


def energy2brho(energy, mass, Q):
    return beta2brho(energy2beta(energy), mass, Q)


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


def input_projectile(comment="", initZ = 70, initA = 175, initAoZ = 2.5, initEnergy = 300.0, initBrho = 6.7, use_url_params = False, initChargeState = 0):
    if use_url_params:
        def _query_value(key):
            value = st.query_params.get(key)
            if isinstance(value, list):
                return value[0] if len(value) > 0 else None
            return value

        def _bounded_int(value, default, min_value, max_value):
            try:
                return max(min_value, min(max_value, int(value)))
            except (TypeError, ValueError):
                return default

        def _bounded_float(value, default, min_value, max_value):
            try:
                return max(min_value, min(max_value, float(value)))
            except (TypeError, ValueError):
                return default

        initZ = _bounded_int(_query_value("z"), initZ, 30, 94)
        initA = _bounded_int(_query_value("a"), initA, 50, 300)
        initAoZ = _bounded_float(_query_value("aoz"), initAoZ, 1.5, 3.6)
        initEnergy = _bounded_float(_query_value("energy"), initEnergy, 50.0, 1000.0)
        initBrho = _bounded_float(_query_value("brho"), initBrho, 2.0, 20.0)
        initChargeState = _bounded_int(_query_value("charge_state"), initChargeState, 0, 6)
        initChargeState = _bounded_int(_query_value("dq"), initChargeState, 0, 6)

    st.write("**Projectile**: Configure the projectile parameters as the initial state for the calculations.")
    if comment != "":
        st.write(comment)

    col1, col2, col3, col4, col5 = st.columns([1, 1, 1, 1, 1.2])
    with col2:
        projectile_Z = st.number_input("Atomic Number (Z)", min_value=30, max_value=94, step=1, value=initZ)
        st.write(f"Element: {z2symbol[projectile_Z]}")
    with col1:
        A = st.number_input("Mass Number (A)", min_value=50, max_value=300, step=1, value=initA)
        st.write(f"A/Z={A/projectile_Z:.4f}")
    with col3:
        charge_states = {0: "Full-strip", 1: "H-like", 2: "He-like", 3: "Li-like", 4: "Be-like", 5: "B-like", 6: "C-like"}
        charge_state = st.selectbox("Charge states", options=list(charge_states.keys()), index=initChargeState, format_func=lambda x: charge_states[x])
        Q = projectile_Z-charge_state
        mass = get_ion_mass_ame20(A, projectile_Z, Q)
        st.write(f"Q={Q}+")
        st.write(f"A/Q={A/Q:.4f}")
    with col4:
        energy_unit = st.selectbox("Energy unit", ["MeV/u", "Tm", "Beta"])
        st.write(f"{mass:.5f} amu")
    with col5:
        if energy_unit == "MeV/u":
            energy = st.number_input("Energy (MeV/u)", min_value=50.0, max_value=1000.0, step=5.0, value=initEnergy)
        elif energy_unit == "Tm":
            brho = st.number_input("Brho (Tm)", min_value=2.0, max_value=20.0, step=0.1, value=initBrho, format="%.3f")
        elif energy_unit == "Beta":
            beta = st.number_input("Beta", min_value=0.3, max_value=0.9, step=0.02, value=0.7, format="%.3f")
        if energy_unit == "MeV/u":
            st.write(f"{energy2brho(energy, mass, Q):.5f} Tm")
            st.write(f"Beta: {energy2beta(energy):.7f}")
        elif energy_unit == "Tm":
            energy = brho2energy(brho, mass, Q)
            st.write(f"{energy:.4f} MeV/u")
            st.write(f"Beta: {energy2beta(energy):.7f}")
        elif energy_unit == "Beta":
            energy = beta2energy(beta)
            st.write(f"{energy:.4f} MeV/u")
            st.write(f"{energy2brho(energy, mass, Q):.5f} Tm")


    return projectile_Z, energy, A, mass, charge_state


def input_materials():
    st.write("**Material Setting**: Select and configure materials. Press 'Add' to append them to the material list.")


    # 選択された物質のリスト (セッションステートを使用)
    if "selected_materials" not in st.session_state:
        st.session_state.selected_materials = ["P10 gas IC-0"]
    if "j" not in st.session_state:
        st.session_state.j = 1

    # 横に並べる
    col1, col2, col3, col4, col5 = st.columns([3, 2, 1.5, 2, 1])

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
        material_units = ["mm", "μm", "cm", "mg/cm2"]
        material_unit = st.radio("Unit", material_units, key=f"material_unit")

    with col4:
        if category == "Gas detectors":
            thickness = st.number_input("Thickness (N/A)", disabled=True)  # 入力不可
            expanded_materials = get_expanded_materials([material])
            st.write(", ".join(expanded_materials))
        elif material_unit == "mg/cm2":
            thickness_mgcm3 = st.number_input("Thickness (mg/cm2)", min_value=0.0, max_value=10000.0, step=1.0, value=100.0)
            st.session_state.thickness = thickness_mgcm3 / density / 1000 / 0.1
        else:
            if material_unit == "μm":
                thickness = st.number_input("Thickness (µm)", min_value=0.0, max_value=10000.0, step=10.0, value=10.0)
                st.session_state.thickness = thickness * 0.001
            elif material_unit == "mm":
                thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=10000.0, step=1.0, value=1.0)
                st.session_state.thickness = thickness * 1
            elif material_unit == "cm":
                thickness = st.number_input("Thickness (cm)", min_value=0.1, max_value=10000.0, step=1.0, value=1.0)
                st.session_state.thickness = thickness * 10
            if category == "Gas":
                st.write(f"{density*1000 * st.session_state.thickness*0.1:.6g} mg/cm2")
            else:
                st.write(f"{density * st.session_state.thickness*0.1:.6g} g/cm2")

    with col5:
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


@lru_cache(maxsize=1)
def load_ame20():
    ame20_data = {}
    import os

    # Download from https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
    file_path = os.path.join(os.path.dirname(__file__), "mass_1.mas20.txt")
    with open(file_path, "r") as f:
        print(file_path, "is loading...")
        lines = f.readlines()[36:]
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        try:
            assert line[4] == " "
            assert line[9] == " "
            assert line[14] == " "
            assert line[27] == " "
            assert line[42] == " "
            Z = int(line[10:14])
            N = int(line[5:9])
            A = Z + N
            mass_excess_keV = float(line[28:42].strip())
            ame20_data[(A, Z)] = mass_excess_keV / 1000.0  # keV -> MeV
        except (ValueError, AssertionError):
            continue
    return ame20_data


mH_MeV = 938.78307389  # hydrogen atom mass energy, MeV
mn_MeV = 939.56542052  # neutron mass energy, MeV


def pairing_delta(A, Z):
    N = A - Z
    if A % 2 == 1:
        return 0.0
    if Z % 2 == 0 and N % 2 == 0:
        return +12.0 / A**0.5
    return -12.0 / A**0.5


def binding_energy_semi_empirical(A, Z):
    av, aS, aC, aA = 15.75, 17.8, 0.711, 23.7
    return av * A - aS * A ** (2.0 / 3.0) - aC * Z * (Z - 1) / A ** (1.0 / 3.0) - aA * (A - 2 * Z) ** 2 / A + pairing_delta(A, Z)


def mass_excess_semi_empirical_mev(A, Z):
    N = A - Z
    B = binding_energy_semi_empirical(A, Z)
    atomic_mass_energy = Z * mH_MeV + N * mn_MeV - B
    mass_excess = atomic_mass_energy - A * mnucleon
    return mass_excess


@lru_cache(maxsize=1)
def load_binding_energy():
    import csv
    import os

    file_path = os.path.join(os.path.dirname(__file__), "binding_energy.csv")
    binding_energy = {}
    with open(file_path, "r", encoding="utf-8", newline="") as f:
        print(file_path, "is loading...")
        for row in csv.DictReader(f):
            Z = int(row["Z"])
            Q = int(row["Q"])
            binding_energy[(Z, Q)] = float(row["Binding Energy (eV)"]) * 1e-6
    return binding_energy


def get_binding_energy(Z, Q):
    if Q == Z:
        return 0.0
    binding_energy = load_binding_energy()
    key = (Z, Q)
    if key not in binding_energy:
        raise ValueError(f"Z={Z}, Q={Q} not found.")
    return binding_energy[key]


def get_ion_mass_ame20(A, Z, Q):
    N = A - Z
    if A <= 0 or Z <= 0 or N < 0 or Q < 0 or Q > Z:
        raise ValueError(f"Invalid nucleus: A={A}, Z={Z}, Q={Q}")

    ame20_data = load_ame20()
    key = (A, Z)
    if key in ame20_data:
        mass_excess = ame20_data[key]
    else:
        mass_excess = mass_excess_semi_empirical_mev(A, Z)
    mass_amu = A + (mass_excess / mnucleon)
    binding_atom = get_binding_energy(Z, 0)
    binding_ion = get_binding_energy(Z, Q)
    binding_correction_amu = (binding_atom - binding_ion) / mnucleon
    return mass_amu - Q * emass / mnucleon + binding_correction_amu
