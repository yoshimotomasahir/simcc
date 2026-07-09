import sys

import pandas as pd
import streamlit as st

sys.path.append("..")
from transmittance_calculation import calculate_transmittance, create_path_table
from utils import *

st.set_page_config(page_title="Transmittance -SimCC-", page_icon="🌠")
st.header("Transmittance due to charge state changes")
st.write("Calculate transmittance through user-defined materials at F0, F1, F3, and F5.")

preset_name = st.sidebar.selectbox("Configuration preset", list(configuration_presets.keys()), key="transmittance_preset")
preset = configuration_presets[preset_name]


if st.session_state.get("transmittance_loaded_preset") != preset_name:
    preset_projectile = preset["projectile"]
    st.session_state["transmittance_projectile_Z"] = preset_projectile["Z"]
    st.session_state["transmittance_projectile_A"] = preset_projectile["A"]
    st.session_state["transmittance_projectile_charge_state"] = preset_projectile["charge_state"]
    st.session_state["transmittance_projectile_energy_unit"] = "Tm"
    st.session_state["transmittance_projectile_brho"] = preset_projectile["brho"]
    st.session_state.pop("transmittance_selected_materials", None)
    st.session_state.pop("transmittance_j", None)
    st.session_state["transmittance_loaded_preset"] = preset_name

projectile = input_projectile(
    use_url_params=True,
    projectile=preset["projectile"],
    key_prefix="transmittance_projectile",
)
projectile_Z = projectile["Z"]
energyD1 = projectile["energy"]
A = projectile["A"]
mass = projectile["mass"]
charge_state = projectile["charge_state"]

exp_correction = input_exp_correction()

materials_by_stage = input_materials(
    stages_with_materials=preset["materials"],
    key_prefix="transmittance",
)

transmittance = calculate_transmittance(projectile, materials_by_stage, exp_correction)
energies = transmittance["energies"]

st.write("**Energy** [MeV/u]")
energy_df = pd.DataFrame(
    [[energies["Initial"], energies["D1"], energies["D2"], energies["D3-D4"], energies["D5-D6"]]],
    columns=["Initial", "D1", "D2", "D3-D4", "D5-D6"],
)
st.dataframe(energy_df.round(4), width="stretch")

st.write("**Transmittance**")
st.dataframe(create_path_table(transmittance["paths"]), width="stretch")

st.write("**Brho** Tm")
brho_rows = []
brho_energies = [energies["D1"], energies["D2"], energies["D3-D4"], energies["D5-D6"]]
for dq in range(7):
    q = projectile_Z - dq
    row = []
    for energy in brho_energies:
        ion_mass = get_ion_mass_ame20(A, projectile_Z, q)
        row.append(energy2brho(energy, ion_mass, q))
    brho_rows.append(row)
brho_df = pd.DataFrame(
    brho_rows,
    index=[f"Z-q={dq}" for dq in range(7)],
    columns=["D1", "D2", "D3-D4", "D5-D6"],
)
st.dataframe(brho_df.round(5), width="stretch")
