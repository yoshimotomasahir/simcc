import streamlit as st

import sys

sys.path.append("..")
from utils import *

st.set_page_config(page_title="Transmittance -SimCC-", page_icon="🌠")
st.header("Transmittance due to charge state changes")
st.write("Calculate transmittance from F0 to the F7 entrance in BigRIPS separator.")

projectile_Z, energy, A, charge_state = input_projectile("Enter the charge state and energy ​​in **D1**.")
dQD1 = charge_state
QD1 = projectile_Z - dQD1
brhoD1 = energy2brho(energy, A, projectile_Z - charge_state)

exp_correction=input_exp_correction()

st.write("**Beamline components**")

col1, col2, col3, col4 = st.columns(4)

with col1:
    target_material = st.text_input("Target element", value="Be")
    target_thickness = st.number_input("Target thickness [mm]", value=1.0, format="%.3f", step=1.0, min_value=0.0)

with col2:
    degrader1_material = st.text_input("DegF1 element", value="Al")
    degrader1_thickness = st.number_input("DegF1 thickness [mm]", value=1.0, format="%.3f", step=1.0, min_value=0.0)

with col3:
    pla3_thickness = st.number_input("Pla3 thickness [mm]", value=0.1, format="%.3f", step=0.1, min_value=0.1)
    st.write("PPAC x2")

with col4:
    degrader2_material = st.text_input("DegF5 element", value="Al")
    degrader2_thickness = st.number_input("DegF5 thickness[mm]", value=1.0, format="%.3f", step=1.0, min_value=0.0)
    degrader2_position = st.radio("DegF5 ** PPAC2", ["before", "after"])
    st.write("PPAC x2")

calculation_option = st.radio("Calculation mode:", ["Calculate charge state transience (excluding F0)", "Assume equilibrium at the last material"])

if calculation_option.startswith("Calculate"):
    calculation_option = "transience"
elif calculation_option.startswith("Assume"):
    calculation_option = "equilibrium"
else:
    raise ValueError("Unknown calculation option selected.")

# D1 energy
energyD1 = brho2energy(brhoD1, A, QD1)


def find_F0_eloss(A, projectile_Z, energyD1, target_material, target_thickness, elossF0a):
    elossF0b = elossF0a
    tolerance = 1e-6
    max_iterations = 1000
    for _ in range(max_iterations):
        elossF0c = GetAnalyticalEloss(A, projectile_Z, energyD1 + elossF0b, target_material, target_thickness)[0]
        if abs(elossF0c - elossF0b) < tolerance:
            return elossF0c
        elossF0b = elossF0c
    return elossF0c


# F0 charge state distribution
elossF0a = GetAnalyticalEloss(A, projectile_Z, energyD1, target_material, target_thickness * 0.1)[0]
elossF0 = find_F0_eloss(A, projectile_Z, energyD1, target_material, target_thickness * 0.1, elossF0a)
MFPF0_entrance = GetMFP(zp=projectile_Z, energy=energyD1 + elossF0, material=target_material, exp_correction=exp_correction)
MFPF0_center = GetMFP(zp=projectile_Z, energy=energyD1 + elossF0 / 2, material=target_material, exp_correction=exp_correction)
P0F0_initial = GetAnalyticalEqProb(MFPF0_entrance)
P0F0 = GetAnalyticalProb(MFPF0_center, target_thickness * 0.1, charge_state=P0F0_initial)

# F1 charge state distribution
elossF1 = GetAnalyticalEloss(A, projectile_Z, energyD1, degrader1_material, degrader1_thickness * 0.1)[0]
MFPF1 = GetMFP(zp=projectile_Z, energy=energyD1 - elossF1 / 2, material=degrader1_material, exp_correction=exp_correction)
P0F1s = []
for dQ in range(7):
    if calculation_option == "transience":
        P0 = np.zeros(7)
        P0[dQ] = 1
        P0 = GetAnalyticalProb(MFPF1, degrader1_thickness * 0.1, charge_state=P0)
    elif calculation_option == "equilibrium":
        P0 = GetAnalyticalEqProb(MFPF1)
    P0F1s.append(P0)
energyD2 = energyD1 - elossF1

# F3 charge state distribution
elossF3Pla = GetAnalyticalEloss(A, projectile_Z, energyD2, "Pla", pla3_thickness * 0.1)[0]
elossF3PPAC1 = GetAnalyticalEloss(A, projectile_Z, energyD2 - elossF3Pla, "Mylar", 0.048 * 0.1)[0]
elossF3PPAC2 = GetAnalyticalEloss(A, projectile_Z, energyD2 - elossF3Pla - elossF3PPAC1, "Mylar", 0.048 * 0.1)[0]
MFPF3Pla = GetMFP(zp=projectile_Z, energy=energyD2 - elossF3Pla / 2, material="Pla", exp_correction=exp_correction)
MFPF3PPAC1 = GetMFP(zp=projectile_Z, energy=energyD2 - elossF3Pla - elossF3PPAC1 / 2, material="Mylar", exp_correction=exp_correction)
MFPF3PPAC2 = GetMFP(zp=projectile_Z, energy=energyD2 - elossF3Pla - elossF3PPAC1 - elossF3PPAC2 / 2, material="Mylar", exp_correction=exp_correction)
P0F3s = []
for dQ in range(7):
    if calculation_option == "transience":
        P0 = np.zeros(7)
        P0[dQ] = 1
        P0 = GetAnalyticalProb(MFPF3Pla, pla3_thickness * 0.1, charge_state=P0)
        GetAnalyticalProb(MFPF3PPAC1, 0.048 * 0.1, charge_state=P0)
        GetAnalyticalProb(MFPF3PPAC2, 0.048 * 0.1, charge_state=P0)
    elif calculation_option == "equilibrium":
        P0 = GetAnalyticalEqProb(MFPF3PPAC2)
    P0F3s.append(P0)
energyD34 = energyD2 - elossF3Pla - elossF3PPAC1 - elossF3PPAC2

# F5 charge state distribution
elossF5PPAC1 = GetAnalyticalEloss(A, projectile_Z, energyD34, "Mylar", 0.048 * 0.1)[0]
MFPF5PPAC1 = GetMFP(zp=projectile_Z, energy=energyD34 - elossF5PPAC1 / 2, material="Mylar", exp_correction=exp_correction)
if degrader2_position == "after":
    elossF5PPAC2 = GetAnalyticalEloss(A, projectile_Z, energyD34 - elossF5PPAC1, "Mylar", 0.048 * 0.1)[0]
    MFPF5PPAC2 = GetMFP(zp=projectile_Z, energy=energyD34 - elossF5PPAC1 - elossF5PPAC2 / 2, material="Mylar", exp_correction=exp_correction)
    elossF5Deg = GetAnalyticalEloss(A, projectile_Z, energyD34 - elossF5PPAC1 - elossF5PPAC2, degrader2_material, degrader2_thickness * 0.1)[0]
    MFPF5Deg = GetMFP(zp=projectile_Z, energy=energyD34 - elossF5PPAC1 - elossF5PPAC2 - elossF5Deg / 2, material=degrader2_material, exp_correction=exp_correction)
else:
    elossF5Deg = GetAnalyticalEloss(A, projectile_Z, energyD34 - elossF5PPAC1, degrader2_material, degrader2_thickness * 0.1)[0]
    MFPF5Deg = GetMFP(zp=projectile_Z, energy=energyD34 - elossF5PPAC1 - elossF5Deg / 2, material=degrader2_material, exp_correction=exp_correction)
    elossF5PPAC2 = GetAnalyticalEloss(A, projectile_Z, energyD34 - elossF5PPAC1 - elossF5Deg, "Mylar", 0.048 * 0.1)[0]
    MFPF5PPAC2 = GetMFP(zp=projectile_Z, energy=energyD34 - elossF5PPAC1 - elossF5Deg - elossF5PPAC2 / 2, material="Mylar", exp_correction=exp_correction)
P0F5s = []
for dQ in range(7):
    if calculation_option == "transience":
        P0 = np.zeros(7)
        P0[dQ] = 1
        if degrader2_position == "after":
            P0 = GetAnalyticalProb(MFPF5PPAC1, 0.048 * 0.1, charge_state=P0)
            P0 = GetAnalyticalProb(MFPF5PPAC2, 0.048 * 0.1, charge_state=P0)
            P0 = GetAnalyticalProb(MFPF5Deg, degrader2_thickness * 0.1, charge_state=P0)
        else:
            P0 = GetAnalyticalProb(MFPF5PPAC1, 0.048 * 0.1, charge_state=P0)
            P0 = GetAnalyticalProb(MFPF5Deg, degrader2_thickness * 0.1, charge_state=P0)
            P0 = GetAnalyticalProb(MFPF5PPAC2, 0.048 * 0.1, charge_state=P0)
    elif calculation_option == "equilibrium":
        if degrader2_position == "after":
            P0 = GetAnalyticalEqProb(MFPF5Deg)
        else:
            P0 = GetAnalyticalEqProb(MFPF5PPAC2)
    P0F5s.append(P0)
energyD56 = energyD34 - elossF5Deg - elossF5PPAC1 - elossF5PPAC2


# Calculate combination of charge states
def compute_path_P0(P0F0, P0F1s, P0F3s, P0F5s):
    P0F0 = np.array(P0F0)
    P0F1s = np.array(P0F1s)
    P0F3s = np.array(P0F3s)
    P0F5s = np.array(P0F5s)

    paths = []
    for i in range(7):
        for j in range(7):
            for k in range(7):
                for l in range(7):
                    prob = P0F0[i] * P0F1s[i, j] * P0F3s[j, k] * P0F5s[k, l]
                    paths.append(((i, j, k, l), prob, [P0F0[i], P0F1s[i, j], P0F3s[j, k], P0F5s[k, l]]))

    paths.sort(key=lambda x: x[1], reverse=True)
    return paths


import pandas as pd

columns = ["D1 [Tm]", "D2 [Tm]", "D34 [Tm]", "D56 [Tm]"]

# Energy
data = [[eval(f"energy{col.split()[0]}") for col in columns]]
df = pd.DataFrame(data, columns=[col.replace("Tm", "MeV/u") for col in columns])
st.dataframe(df, width="stretch")

# Transmittance
dQ2label = {
    0: "Full",
    1: "H-",
    2: "He-",
    3: "Li-",
    4: "Be-",
    5: "B-",
    6: "C-",
}


def create_path_table(path_P0):
    data = []
    for path, prob, probs in path_P0:
        if prob < 0.001:
            break
        row = {
            "Total": f"{prob:.2%}",
            "Z-q=": ",".join(map(str, path)),
            "D1": dQ2label[path[0]],
            "p1": f"{probs[0]:.2%}",
            "D2": dQ2label[path[1]],
            "p2": f"{probs[1]:.2%}",
            "D3D4": dQ2label[path[2]],
            "p3": f"{probs[2]:.2%}",
            "D5D6": dQ2label[path[3]],
            "p4": f"{probs[3]:.2%}",
        }
        data.append(row)

    df = pd.DataFrame(data)
    return df


df_paths = create_path_table(compute_path_P0(P0F0, P0F1s, P0F3s, P0F5s))
st.dataframe(df_paths, width="stretch")

# Brho
data = []
columns = ["D1 [Tm]", "D2 [Tm]", "D34 [Tm]", "D56 [Tm]"]
for dQ in range(7):
    row = [eval(f"energy2brho(energy{col.split()[0]}, A, projectile_Z-dQ)") for col in columns]
    data.append(row)
index = [f"Z-q={dQ}" for dQ in range(7)]
df = pd.DataFrame(data, index=index, columns=columns)
st.dataframe(df, width="stretch")
