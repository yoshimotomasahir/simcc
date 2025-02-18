import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
from simcc import *
import sys
sys.path.append("..")
from utils import *

st.title("Energy resolution")

st.write("Projectile:")
col1, col2, col3 = st.columns(3)
with col1:
    A = st.number_input("Mass number A", min_value=1, max_value=300, value=238, step=1)
with col2:
    Z = st.number_input("Atomic number Z", min_value=1, max_value=120, value=92, step=1)
with col3:
    Q = st.number_input("Charge Q", min_value=1, max_value=120, value=90, step=1)

col4, col5, col6 = st.columns(3)
with col4:
    energy0 = st.number_input("Min Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=150.0, step=5.0)
with col5:
    energy1 = st.number_input("Max Energy [MeV/u]", min_value=100.0, max_value=1000.0, value=350.0, step=5.0)
with col6:
    energyBin = st.number_input("Number of bins", min_value=2, max_value=100, value=5, step=1)

st.write("Materials:")
materials = []
col1, col2, col3 = st.columns(3)
i = 0
with col1:
    material = st.selectbox(
        f"Material",
        options=materialOptions,
        index=5,
        key=f"material_{i}",
    )
with col2:
    thickness = st.number_input(
        "Thickness",
        min_value=0.0,
        max_value=10000.0,
        value=480.0,
        step=1.0,
        key=f"thickness_{i}",
        format="%.3f",
    )
with col3:
    thicknessUnit = st.selectbox(
        f"Thickness unit",
        options=["mm", "mg/cm2"],
        index=0,
        key=f"charge_state_{i}",
    )
materials.append({"Material": material, "Thickness": thickness, "ThicknessUnit": thicknessUnit})

if st.button("Execute Calculation"):
    st.write("### Calculation Results")

    Eres = {}
    Zres = {}
    energies = np.linspace(energy0, energy1, energyBin)
    for energy in energies:
        Ein = energy
        for i, material in enumerate(materials, start=1):
            if material["Thickness"] == 0:
                continue

            density = GetMaterial(material["Material"])["density"]
            if material["ThicknessUnit"] == "mg/cm2":
                thickness = material["Thickness"] / density * 0.010
            thicknessStr = f"{thickness*1000:.2f} μm" if thickness < 0.1 else f"{thickness:.2f} mm"

            st.write(f"#### $^{{{A}}}${z2symbol[Z]}$^{{{Q}+}}$ {Ein:.1f} MeV/u into {material['Material']} {thicknessStr}")
            st.write(f"density: " + (f"{density:.2f} g/cm3" if density > 0.1 else f"{density*1000:.4f} mg/cm3"))

            histories = None
            dEtotal, dEcol, dEcc, charges, histories = GetMCEloss(
                A,
                Z,
                Q,
                Ein,
                material["Material"],
                thickness * 0.1,
                histories=histories,
            )
            histories = None
            dEtotalm1, dEcolm1, dEccm1, chargesm1, historiesm1 = GetMCEloss(
                A,
                Z - 1,
                Q - 1,
                Ein,
                material["Material"],
                thickness * 0.1,
                histories=histories,
            )

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
            for param, dE, dEm1 in zip(["dEcol", "dEtotal (col+cc)"], [dEcol, dEtotal], [dEcolm1, dEtotalm1]):
                label = f"{param}\nMean:{np.mean(dE):.1f}\nStdev:{np.std(dE):.2f} ({np.std(dE)/np.mean(dE):.2%})"
                labelm1 = f"{param}(Z-1)\nMean:{np.mean(dEm1):.1f}\nStdev:{np.std(dEm1):.2f} ({np.std(dEm1)/np.mean(dEm1):.2%})"
                ax.hist(dE, alpha=0.5, bins=50, range=histRange, label=label)
                if param not in Eres:
                    Eres[param] = []
                    Zres[param] = []
                Eres[param].append(np.std(dE) / np.mean(dE))
                Zres[param].append(np.std(dE) / (np.mean(dE) - np.mean(dEm1)))
                if "total" in label:
                    ax.hist(dEm1, alpha=0.5, bins=50, range=histRange, label=labelm1)
            ax.legend(loc="upper left", bbox_to_anchor=(0.7, 1))
            ax.set_xlabel("Energy loss [MeV/u]")
            ax.set_title("Energy loss distribution")
            ax.set_xlim(*histRange)
            st.pyplot(fig)

    fig, ax = plt.subplots(figsize=(5, 2))
    for key, val in Eres.items():
        ax.plot(energies, np.array(val) * 100, label=key)
    ax.set_xlabel("Energy [MeV/u]")
    ax.set_ylabel("Energy resolution (σ) [%]")
    ax.set_ylim(0, 1.6)
    ax.grid(alpha=0.4)
    ax.legend()
    st.pyplot(fig)

    fig, ax = plt.subplots(figsize=(5, 2))
    for key, val in Zres.items():
        ax.plot(energies, 1 / np.array(val), label=key)
    ax.set_xlabel("Energy [MeV/u]")
    ax.set_ylabel("1 / Z resolution (σ)")
    ax.set_ylim(0, 10)
    ax.grid(alpha=0.4)
    ax.legend()
    st.pyplot(fig)

    fig, ax = plt.subplots(figsize=(5, 2))
    for key, val in Zres.items():
        ax.plot(energies, np.array(val), label=key)
    ax.set_xlabel("Energy [MeV/u]")
    ax.set_ylabel("Z resolution (σ)")
    ax.set_ylim(0, 0.5)
    ax.grid(alpha=0.4)
    ax.legend()
    st.pyplot(fig)

    st.success("Calculation executed successfully!")
