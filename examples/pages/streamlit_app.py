# "streamlit run streamlit_app.py" to run
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
from simcc import *

st.write("Projectile:")
col1, col2, col3, col4 = st.columns(4)
with col1:
    A = st.number_input("Mass number A", min_value=1, max_value=300, value=170, step=1)
with col2:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=63, step=1)
with col3:
    Q = st.number_input("Charge Q", min_value=1, max_value=120, value=63, step=1)
with col4:
    energy = st.number_input(
        "Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=300.0, step=5.0
    )

st.write("Materials:")
materials = []
initialMaterial = {1: 0, 2: 1, 3: 2, 4: 3, 5: 2, 6: 2}
initialThickness = {1: 2.0, 2: 1.0, 3: 0.048, 4: 0.1, 5: 0.048, 6: 0.096}
initialChargeState = {1: 0, 2: 2, 3: 0, 4: 0, 5: 0, 6: 0}
for i in range(1, 7):
    col1, col2, col3 = st.columns(3)
    with col1:
        material = st.selectbox(
            f"Material",
            options=[
                "Be",
                "Al",
                "Mylar",
                "Pla",
                "Kapton",
                "P10",
                "Xe7",
                "Diamond",
                "Gold",
            ],
            index=initialMaterial[i],
            key=f"material_{i}",
        )
    with col2:
        thickness = st.number_input(
            "Thickness [mm]",
            min_value=0.0,
            max_value=1000.0,
            value=initialThickness[i],
            step=1.0,
            key=f"thickness_{i}",
            format="%.3f",
        )
    with col3:
        chargeState = st.selectbox(
            f"Charge state selected after passage",
            options=[
                "All",
                "Full-strip",
                "H-like",
                "He-like",
            ],
            index=initialChargeState[i],
            key=f"charge_state_{i}",
        )
    materials.append(
        {"Material": material, "Thickness": thickness, "ChargeState": chargeState}
    )

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    Ein = energy
    Transmission = 1.0
    for i, material in enumerate(materials, start=1):
        if material["Thickness"] == 0:
            continue

        thicknessStr = (
            f"{material['Thickness']*1000} μm"
            if material["Thickness"] < 0.1
            else f"{material['Thickness']} mm"
        )

        st.write(
            f"#### $^{{{A}}}${z2symbol[Z]}$^{{{Q}+}}$ {Ein:.1f} MeV/u into {material['Material']} {thicknessStr}"
        )

        if Q != "All":
            histories = None

        dEtotal, dEcol, dEcc, charges, histories = GetDeltaE(
            A,
            Z,
            Q,
            Ein,
            material["Material"],
            material["Thickness"] * 0.1,
            histories=histories,
        )
        Ein = Ein - np.mean(dEtotal)
        fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
        ax = axes[0]
        for q in range(Z, Z - 4, -1):
            ax.plot(charges["length"] * 10, charges[q], label=f"Q={q}")
        ax.legend()
        ax.set_xlabel("Thickness [mm]")
        ax.set_xlim(0, None)
        ax.set_ylim(0, 1)
        ax.set_title("Charge state fraction")

        ax = axes[1]
        histRange = [np.mean(dEtotal) * 0.95, np.mean(dEtotal) * 1.05]
        for param, dE in zip(["dEcol", "dEtotal (col+cc)"], [dEcol, dEtotal]):
            label = f"{param}\nMean:{np.mean(dE):.1f}\nStdev:{np.std(dE):.2f} ({np.std(dE)/np.mean(dE):.2%})"
            ax.hist(dE, alpha=0.5, bins=50, range=histRange, label=label)
        ax.legend(loc="upper left", bbox_to_anchor=(0.7, 1))
        ax.set_xlabel("Energy loss [MeV/u]")
        ax.set_title("Energy loss distribution")
        ax.set_xlim(*histRange)
        st.pyplot(fig)

        if material["ChargeState"] == "All":
            Q = material["ChargeState"]
        elif material["ChargeState"] == "Full-strip":
            Q = Z
            Transmission *= charges[Q][-1]
        elif material["ChargeState"] == "H-like":
            Q = Z - 1
            Transmission *= charges[Q][-1]
        elif material["ChargeState"] == "He-like":
            Q = Z - 2
            Transmission *= charges[Q][-1]
        st.write(f"Transmission: {Transmission:.2%}")

    st.success("Calculation executed successfully!")
