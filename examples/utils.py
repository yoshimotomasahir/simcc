import streamlit as st
import numpy as np

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
    "Strippers": ["Ta", "Al", "W", "Pt", "Au"],
    "Gas detectors": ["P10 gas IC", "Xe gas IC", "PPAC"],
    "Gas": ["P10", "Xe7", "iC4H10"],
    "Other detectors": ["Plastic", "Diamond"],
    "Degraders": ["Al", "Cu"],
    "Targets": ["Be", "W"],
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
    st.write("**Projectile**")

    col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
    with col1:
        projectile_Z = st.slider("Projectile Z", 30, 92, 70)
    with col2:
        energy = st.slider("Energy [MeV/u]", 50, 1000, 250, step=5)
    with col3:
        AoZ = st.slider("A/Z", 1.5, 3.5, 2.5, step=0.01)
        A = int(np.round(projectile_Z * AoZ))
        st.write(f"A={A}")
    with col4:
        charge_states = {0: "Full-strip", 1: "H-like", 2: "He-like", 3: "Li-like", 4: "Be-like", 5: "B-like", 6: "C-like"}
        charge_state = st.selectbox("Charge states", options=list(charge_states.keys()), format_func=lambda x: charge_states[x])
        st.write(f"Brho={energy2brho(energy, A, projectile_Z-charge_state):.5f} Tm")

    return projectile_Z, energy, A, charge_state


def input_materials():
    st.write("**Material**")

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
        elif category == "Gas":
            thickness = st.number_input("Thickness (mm)", min_value=0.1, max_value=5000.0, step=0.1, value=600.0)
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
    return st.session_state.selected_materials
