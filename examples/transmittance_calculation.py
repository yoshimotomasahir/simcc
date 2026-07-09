import numpy as np
import pandas as pd

from simcc import GetAnalyticalEloss, GetAnalyticalEqProb, GetAnalyticalProb, GetMFP
from utils import get_expanded_materials, get_material_name_length


def expand_stage_materials(materials_by_stage, stage):
    return get_expanded_materials(materials_by_stage.get(stage, []))


def last_material_name(materials):
    expanded_materials = get_expanded_materials(materials)
    if len(expanded_materials) == 0:
        return None
    return get_material_name_length(expanded_materials[-1])[0]


def equilibrium_probability(projectile_z, charge_state, material_name, energy, exp_correction):
    if material_name is None:
        p0 = np.zeros(7)
        p0[charge_state] = 1
        return p0
    mfp = GetMFP(zp=projectile_z, energy=energy, material=material_name, exp_correction=exp_correction)
    return GetAnalyticalEqProb(mfp)


def stage_energy_out(a, projectile_z, energy_in, materials):
    energy = energy_in
    for material in materials:
        material_name, length = get_material_name_length(material)
        eloss = GetAnalyticalEloss(a, projectile_z, energy, material_name, length * 0.1, z_effective=1)[0]
        energy -= eloss
    return energy


def find_energy_in_for_material(a, projectile_z, energy_out, material_name, length):
    eloss = GetAnalyticalEloss(a, projectile_z, energy_out, material_name, length * 0.1)[0]
    tolerance = 1e-6
    max_iterations = 1000
    for _ in range(max_iterations):
        next_eloss = GetAnalyticalEloss(a, projectile_z, energy_out + eloss, material_name, length * 0.1)[0]
        if abs(next_eloss - eloss) < tolerance:
            return energy_out + next_eloss
        eloss = next_eloss
    return energy_out + eloss


def stage_energy_in(a, projectile_z, energy_out, materials):
    energy = energy_out
    for material in reversed(materials):
        material_name, length = get_material_name_length(material)
        energy = find_energy_in_for_material(a, projectile_z, energy, material_name, length)
    return energy


def stage_transition_matrix(a, projectile_z, energy_in, materials, exp_correction):
    matrix = []
    energy_out = stage_energy_out(a, projectile_z, energy_in, materials)

    if len(materials) == 0:
        return np.eye(7), energy_out

    for dq in range(7):
        energy = energy_in
        p0 = np.zeros(7)
        p0[dq] = 1
        for material in materials:
            material_name, length = get_material_name_length(material)
            eloss = GetAnalyticalEloss(a, projectile_z, energy, material_name, length * 0.1, z_effective=1)[0]
            mfp = GetMFP(zp=projectile_z, energy=energy - eloss / 2, material=material_name, exp_correction=exp_correction)
            p0 = GetAnalyticalProb(mfp, length * 0.1, charge_state=p0)
            energy -= eloss
        matrix.append(p0)
    return np.array(matrix), energy_out


def compute_paths(p_d1, t_d2, t_d34, d_d56):
    paths = []
    for dq1 in range(7):
        for dq2 in range(7):
            for dq34 in range(7):
                for dq56 in range(7):
                    probs = [
                        p_d1[dq1],
                        t_d2[dq1, dq2],
                        t_d34[dq2, dq34],
                        d_d56[dq34, dq56],
                    ]
                    total = np.prod(probs)
                    paths.append(((dq1, dq2, dq34, dq56), total, probs))
    paths.sort(key=lambda x: x[1], reverse=True)
    return paths


def top_paths(paths, n=5):
    return paths[:n]


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


def calculate_transmittance(projectile, materials_by_stage, exp_correction):
    projectile_z = projectile["Z"]
    energy_d1 = projectile["energy"]
    a = projectile["A"]
    charge_state = projectile["charge_state"]

    f0_materials = expand_stage_materials(materials_by_stage, "F0")
    f1_materials = expand_stage_materials(materials_by_stage, "F1")
    f3_materials = expand_stage_materials(materials_by_stage, "F3")
    f5_materials = expand_stage_materials(materials_by_stage, "F5")

    energy_initial = stage_energy_in(a, projectile_z, energy_d1, f0_materials)
    p_d1 = equilibrium_probability(
        projectile_z,
        charge_state,
        last_material_name(materials_by_stage["F0"]),
        energy_d1,
        exp_correction,
    )

    t_d2, energy_d2 = stage_transition_matrix(a, projectile_z, energy_d1, f1_materials, exp_correction)
    t_d34, energy_d34 = stage_transition_matrix(a, projectile_z, energy_d2, f3_materials, exp_correction)
    d_d56, energy_d56 = stage_transition_matrix(a, projectile_z, energy_d34, f5_materials, exp_correction)
    paths = compute_paths(p_d1, t_d2, t_d34, d_d56)

    return {
        "energies": {
            "Initial": energy_initial,
            "D1": energy_d1,
            "D2": energy_d2,
            "D3-D4": energy_d34,
            "D5-D6": energy_d56,
        },
        "paths": paths,
        "stage_materials": {
            "F0": f0_materials,
            "F1": f1_materials,
            "F3": f3_materials,
            "F5": f5_materials,
        },
        "probabilities": {
            "F0": p_d1,
            "F1": t_d2,
            "F3": t_d34,
            "F5": d_d56,
        },
    }
