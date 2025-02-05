import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
from simcc import *

st.write("Projectile:")
col1, col2, col3 = st.columns(3)
with col1:
    A = st.number_input("Mass number A", min_value=1, max_value=300, value=119, step=1)
with col2:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=58, step=1)
with col3:
    Q = st.number_input("Charge Q", min_value=1, max_value=120, value=58, step=1)

col4, col5, col6 = st.columns(3)
with col4:
    energy0 = st.number_input(
        "Min Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=150.0, step=5.0
    )
with col5:
    energy1 = st.number_input(
        "Max Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=340.0, step=5.0
    )
with col6:
    energy_bin = st.number_input(
        "Number of bins", min_value=2, max_value=100, value=5, step=1
    )

st.write("Materials:")
materials = []
col1, col2, col3 = st.columns(3)
i = 0
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
        index=1,
        key=f"material_{i}",
    )
with col2:
    thickness = st.number_input(
        "Thickness [mm]",
        min_value=0.0,
        max_value=1000.0,
        value=2.0,
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
        index=0,
        key=f"charge_state_{i}",
    )
materials.append(
    {"Material": material, "Thickness": thickness, "ChargeState": chargeState}
)

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    stds = {}
    energies = np.linspace(energy0, energy1, energy_bin)
    for energy in energies:
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
                if param not in stds:
                    stds[param] = []
                stds[param].append(np.std(dE)/np.mean(dE))
            ax.legend(loc="upper left", bbox_to_anchor=(0.7, 1))
            ax.set_xlabel("Energy loss [MeV/u]")
            ax.set_title("Energy loss distribution")
            ax.set_xlim(*histRange)
            st.pyplot(fig)
            st.write(f"Transmission: {Transmission:.2%}")

    fig, ax = plt.subplots(figsize=(5, 4))
    
    for key,val in stds.items():
        ax.plot(energies, np.array(val)*100, label=key)
    ax.set_xlabel("Energy [MeV/u]")
    ax.set_ylabel("Energy resolution (σ) [%]")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.4)
    ax.legend()
    st.pyplot(fig)

    st.success("Calculation executed successfully!")
