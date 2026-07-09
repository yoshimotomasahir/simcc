import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.append("..")
from simcc import GetAnalyticalEloss, GetMaterial
from transmittance_calculation import calculate_stage_energies, calculate_transmittance, top_paths
from utils import *

st.set_page_config(page_title="Transmittance Optimisation -SimCC-", page_icon="🌠")
st.header("Transmittance Optimisation")
st.write("Compare how candidate stripper materials change the charge-state transmittance.")

preset_name = st.sidebar.selectbox("Configuration preset", list(configuration_presets.keys()), key="material_optimisation2_preset")
preset = configuration_presets[preset_name]

if st.session_state.get("material_optimisation2_loaded_preset") != preset_name:
    preset_projectile = preset["projectile"]
    st.session_state["material_optimisation2_projectile_Z"] = preset_projectile["Z"]
    st.session_state["material_optimisation2_projectile_A"] = preset_projectile["A"]
    st.session_state["material_optimisation2_projectile_charge_state"] = preset_projectile["charge_state"]
    st.session_state["material_optimisation2_projectile_energy_unit"] = "Tm"
    st.session_state["material_optimisation2_projectile_brho"] = preset_projectile["brho"]
    st.session_state.pop("material_optimisation2_selected_materials", None)
    st.session_state.pop("material_optimisation2_j", None)
    st.session_state["material_optimisation2_loaded_preset"] = preset_name

projectile = input_projectile(
    use_url_params=True,
    projectile=preset["projectile"],
    key_prefix="material_optimisation2_projectile",
)
projectile_z = projectile["Z"]
a = projectile["A"]

exp_correction = input_exp_correction()

materials_by_stage = input_materials(
    stages_with_materials=preset["materials"],
    key_prefix="material_optimisation2",
)

baseline_energies = calculate_stage_energies(projectile, materials_by_stage)

st.write("**Candidate stripper**")
stage = st.radio("Stripper insert position", ["F1", "F3", "F5"], index=1, horizontal=True)
dE = st.number_input(
    "Energy loss in the stripper [MeV/u]",
    value=2.0,
    min_value=0.0,
    format="%.3f",
)


def get_thickness(a, z, energy, dE, material, tmax=200.0):
    if energy - dE <= 0:
        return np.nan

    def f(t):
        return GetAnalyticalEloss(a, z, energy, material, t, z_effective=1)[0] - dE

    t_hi = 1e-4
    while t_hi < tmax and f(t_hi) < 0:
        t_hi *= 1.6
    if t_hi >= tmax:
        return np.nan

    from scipy.optimize import root_scalar

    sol = root_scalar(f, bracket=(0.0, t_hi), method="brentq")
    return sol.root if sol.converged else np.nan


energy_by_stage = {
    "F1": baseline_energies["D2"],
    "F3": baseline_energies["D3-D4"],
    "F5": baseline_energies["D5-D6"],
}
stripper_energy_in = energy_by_stage[stage]

elements = ["Carbon", "Al", "Ti", "Cu", "Nb", "Ta"]
thickness_rows = []

for element in elements:
    thickness_cm = get_thickness(a, projectile_z, stripper_energy_in, dE, element)
    result = GetMaterial(element)
    zts, density = result["zts"], result["density"]
    thickness_rows.append(
        {
            "Stripper": element,
            "Z": zts[0],
            "Thickness [um]": thickness_cm * 10000 if np.isfinite(thickness_cm) else np.nan,
            "x [mg/cm2]": thickness_cm * density * 1000 if np.isfinite(thickness_cm) else np.nan,
            "rho [mg/cm3]": density * 1000,
            "thickness_cm": thickness_cm,
        }
    )

thickness_df = pd.DataFrame(thickness_rows)
st.write("**Equivalent thickness for the energy loss**")
st.dataframe(thickness_df.drop(columns=["thickness_cm"]).round(3), width="stretch")


def create_fixed_path_table(paths, target_paths, energy_d56):
    path_map = {path: (total, probs) for path, total, probs in paths}
    rows = []
    for path in target_paths:
        total, probs = path_map[path]
        rows.append(
            {
                "Total": f"{total:.2%}",
                "Z-q": ",".join(map(str, path)),
                "p_F0": f"{probs[0]:.2%}",
                "p_F1": f"{probs[1]:.2%}",
                "p_F3": f"{probs[2]:.2%}",
                "p_F5": f"{probs[3]:.2%}",
                "D5-D6 [MeV/u]": energy_d56,
            }
        )
    return pd.DataFrame(rows)

if st.button("Calculate transmittance change"):
    baseline = calculate_transmittance(projectile, materials_by_stage, exp_correction)
    baseline_top_paths = top_paths(baseline["paths"], 5)
    target_paths = [path for path, _, _ in baseline_top_paths]

    change_by_path = {path: [] for path, _, _ in baseline_top_paths}
    baseline_top_path_table = create_fixed_path_table(baseline["paths"], target_paths, baseline["energies"]["D5-D6"])
    baseline_top_path_table.insert(0, "Stripper", "None")
    trial_top_path_tables = [baseline_top_path_table]

    for thickness_row in thickness_rows:
        element = thickness_row["Stripper"]
        thickness_cm = thickness_row["thickness_cm"]

        if not np.isfinite(thickness_cm):
            for path in change_by_path:
                change_by_path[path].append(np.nan)
            trial_top_path_tables.append(pd.DataFrame([{"Stripper": element}]))
            continue

        added_material = f"{element} {thickness_cm*10:.6g} mm"
        trial_materials_by_stage = {key: value.copy() for key, value in materials_by_stage.items()}
        trial_materials_by_stage[stage].append(added_material)
        trial = calculate_transmittance(projectile, trial_materials_by_stage, exp_correction)
        trial_totals = {path: total for path, total, _ in trial["paths"]}
        trial_top_path_table = create_fixed_path_table(trial["paths"], target_paths, trial["energies"]["D5-D6"])
        trial_top_path_table.insert(0, "Stripper", element)
        trial_top_path_tables.append(trial_top_path_table)

        for path, baseline_total, _ in baseline_top_paths:
            trial_total = trial_totals.get(path, 0.0)
            change = (trial_total / baseline_total - 1.0) * 100 if baseline_total > 0 else np.nan
            change_by_path[path].append(change)

    fig = go.Figure()
    for path, changes in change_by_path.items():
        fig.add_trace(
            go.Scatter(
                x=elements,
                y=changes,
                mode="lines",
                name=f"Z-q={','.join(map(str, path))}",
            )
        )
    fig.update_layout(
        xaxis_title="Added stripper",
        yaxis_title="Transmission change [%]",
    )
    st.write("**Transmittance change after adding the stripper**")
    st.plotly_chart(fig, width="stretch")

    st.write("**Top 5 charge-state paths**")
    st.dataframe(pd.concat(trial_top_path_tables, ignore_index=True).round(4), width="stretch")
