import sys

import numpy as np
import pandas as pd
import streamlit as st

sys.path.append("..")
from utils import *

st.set_page_config(page_title="Transmittance -SimCC-", page_icon="🌠")
st.header("Transmittance due to charge state changes")
st.write("Calculate transmittance through user-defined materials at F0, F1, F3, and F5.")

samples = {
    "195-Ta(HL) 6.7618 Tm": {
        "projectile": {"Z": 73, "A": 195, "brho": 6.7618, "charge_state": 1},
        "materials": {
            "F0": ["Be 5 mm"],
            "F1": ["Al 1.4 mm"],
            "F3": ["Plastic 0.1 mm", "PPAC", "Diamond 0.2 mm", "PPAC"],
            "F5": ["PPAC", "PPAC", "Al 1 mm"],
        },
    },
    "184-Tm(HL) 6.7618 Tm": {
        "projectile": {"Z": 69, "A": 184, "brho": 6.7618, "charge_state": 1},
        "materials": {
            "F0": ["Be 5 mm"],
            "F1": ["Al 1.4 mm"],
            "F3": ["Plastic 0.1 mm", "PPAC", "Diamond 0.2 mm", "PPAC", "Al 0.2 mm"],
            "F5": ["PPAC", "PPAC", "Al 1 mm"],
        },
    },
    "151-Hf(HL) 5.7146 Tm": {
        "projectile": {"Z": 72, "A": 151, "brho": 5.7146, "charge_state": 1},
        "materials": {
            "F0": ["Be 2 mm"],
            "F1": ["Ta 0.01 mm"],
            "F3": ["Plastic 0.1 mm", "PPAC", "Diamond 0.2 mm", "PPAC"],
            "F5": ["PPAC", "Al 2 mm", "PPAC"],
        },
    },
    "116-La(FS) 5.2088 Tm": {
        "projectile": {"Z": 57, "A": 116, "brho": 5.2088, "charge_state": 0},
        "materials": {
            "F0": ["Be 3 mm"],
            "F1": ["Al 1.4 mm"],
            "F3": ["Plastic 0.1 mm", "PPAC", "Diamond 0.2 mm", "PPAC"],
            "F5": ["PPAC", "Al 1.5 mm", "PPAC"],
        },
    },
}

sample_name = st.sidebar.selectbox("Sample setup", list(samples.keys()), key="transmittance_sample")
sample = samples[sample_name]


if st.session_state.get("transmittance_loaded_sample") != sample_name:
    sample_projectile = sample["projectile"]
    st.session_state["transmittance_projectile_Z"] = sample_projectile["Z"]
    st.session_state["transmittance_projectile_A"] = sample_projectile["A"]
    st.session_state["transmittance_projectile_charge_state"] = sample_projectile["charge_state"]
    st.session_state["transmittance_projectile_energy_unit"] = "Tm"
    st.session_state["transmittance_projectile_brho"] = sample_projectile["brho"]
    st.session_state.pop("transmittance_selected_materials", None)
    st.session_state.pop("transmittance_j", None)
    st.session_state["transmittance_loaded_sample"] = sample_name

projectile = input_projectile(
    use_url_params=True,
    projectile=sample["projectile"],
    key_prefix="transmittance_projectile",
)
projectile_Z = projectile["Z"]
energyD1 = projectile["energy"]
A = projectile["A"]
mass = projectile["mass"]
charge_state = projectile["charge_state"]

exp_correction = input_exp_correction()

materials_by_stage = input_materials(
    stages_with_materials=sample["materials"],
    key_prefix="transmittance",
)

def expanded_stage_materials(stage):
    return get_expanded_materials(materials_by_stage.get(stage, []))


def last_material_name(materials):
    expanded_materials = get_expanded_materials(materials)
    if len(expanded_materials) == 0:
        return None
    return get_material_name_length(expanded_materials[-1])[0]


def equilibrium_probability(material_name, energy):
    if material_name is None:
        P0 = np.zeros(7)
        P0[charge_state] = 1
        return P0
    MFP = GetMFP(zp=projectile_Z, energy=energy, material=material_name, exp_correction=exp_correction)
    return GetAnalyticalEqProb(MFP)


def stage_energy_out(energy_in, materials):
    energy = energy_in
    for material in materials:
        material_name, length = get_material_name_length(material)
        eloss = GetAnalyticalEloss(A, projectile_Z, energy, material_name, length * 0.1, z_effective=1)[0]
        energy -= eloss
    return energy


def find_energy_in_for_material(energy_out, material_name, length):
    eloss = GetAnalyticalEloss(A, projectile_Z, energy_out, material_name, length * 0.1)[0]
    tolerance = 1e-6
    max_iterations = 1000
    for _ in range(max_iterations):
        next_eloss = GetAnalyticalEloss(A, projectile_Z, energy_out + eloss, material_name, length * 0.1)[0]
        if abs(next_eloss - eloss) < tolerance:
            return energy_out + next_eloss
        eloss = next_eloss
    return energy_out + eloss


def stage_energy_in(energy_out, materials):
    energy = energy_out
    for material in reversed(materials):
        material_name, length = get_material_name_length(material)
        energy = find_energy_in_for_material(energy, material_name, length)
    return energy


def stage_transition_matrix(energy_in, materials):
    matrix = []
    energy_out = stage_energy_out(energy_in, materials)

    if len(materials) == 0:
        return np.eye(7), energy_out

    for dq in range(7):
        energy = energy_in
        P0 = np.zeros(7)
        P0[dq] = 1
        for material in materials:
            material_name, length = get_material_name_length(material)
            eloss = GetAnalyticalEloss(A, projectile_Z, energy, material_name, length * 0.1, z_effective=1)[0]
            MFP = GetMFP(zp=projectile_Z, energy=energy - eloss / 2, material=material_name, exp_correction=exp_correction)
            P0 = GetAnalyticalProb(MFP, length * 0.1, charge_state=P0)
            energy -= eloss
        matrix.append(P0)
    return np.array(matrix), energy_out


def compute_paths(P_D1, T_D2, T_D34, D_D56):
    paths = []
    for dq1 in range(7):
        for dq2 in range(7):
            for dq34 in range(7):
                for dq56 in range(7):
                    probs = [
                        P_D1[dq1],
                        T_D2[dq1, dq2],
                        T_D34[dq2, dq34],
                        D_D56[dq34, dq56],
                    ]
                    total = np.prod(probs)
                    paths.append(((dq1, dq2, dq34, dq56), total, probs))
    paths.sort(key=lambda x: x[1], reverse=True)
    return paths


def create_path_table(paths):
    rows = []
    for path, total, probs in paths:
        if total < 0.001:
            break
        row = {
            "Total": f"{total:.2%}",
            "Z-q": ",".join(map(str, path)),
            "p_F0": f"{probs[0]:.2%}",
            "p_F1": f"{probs[1]:.2%}",
            "p_F3": f"{probs[2]:.2%}",
            "p_F5": f"{probs[3]:.2%}",
        }
        rows.append(row)
    return pd.DataFrame(rows)


F0_materials = expanded_stage_materials("F0")
F1_materials = expanded_stage_materials("F1")
F3_materials = expanded_stage_materials("F3")
F5_materials = expanded_stage_materials("F5")

energyInitial = stage_energy_in(energyD1, F0_materials)
P_D1 = equilibrium_probability(last_material_name(materials_by_stage["F0"]), energyD1)

T_D2, energyD2 = stage_transition_matrix(energyD1, F1_materials)
T_D34, energyD34 = stage_transition_matrix(energyD2, F3_materials)
D_D56, energyD56 = stage_transition_matrix(energyD34, F5_materials)

st.write("**Energy** [MeV/u]")
energy_df = pd.DataFrame(
    [[energyInitial, energyD1, energyD2, energyD34, energyD56]],
    columns=["Initial", "D1", "D2", "D3-D4", "D5-D6"],
)
st.dataframe(energy_df.round(4), width="stretch")

st.write("**Transmittance**")
st.dataframe(create_path_table(compute_paths(P_D1, T_D2, T_D34, D_D56)), width="stretch")

st.write("**Brho** Tm")
brho_rows = []
energies = [energyD1, energyD2, energyD34, energyD56]
for dq in range(7):
    q = projectile_Z - dq
    row = []
    for energy in energies:
        ion_mass = get_ion_mass_ame20(A, projectile_Z, q)
        row.append(energy2brho(energy, ion_mass, q))
    brho_rows.append(row)
brho_df = pd.DataFrame(
    brho_rows,
    index=[f"Z-q={dq}" for dq in range(7)],
    columns=["D1", "D2", "D3-D4", "D5-D6"],
)
st.dataframe(brho_df.round(5), width="stretch")
