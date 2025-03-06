import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pycatima as catima
from simcc import *
import sys

sys.path.append("..")
from utils import *

import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.collections as mc
import matplotlib.cm as cm

st.set_page_config(page_title="Energy loss distribution -SimCC-", page_icon="🌠")
st.title("MC for energy loss distribution")

projectile_Z, energy, A, charge_state = input_projectile()

materials = input_materials()

if st.button("Execute Calculation"):

    expanded_materials = get_expanded_materials(materials)
    Ein = energy
    histories = None
    total_length = 0
    rs = np.random.RandomState(1)
    rs_mc = np.random.RandomState(1)
    for i, material in enumerate(expanded_materials):
        material_name, length = get_material_name_length(material)

        dEtotal, dEcol, dEcc, charges, histories = GetMCEloss(
            A,
            projectile_Z,
            projectile_Z - charge_state,
            Ein,
            material_name,
            length * 0.1,
            random_state=rs_mc,
            histories=histories,
        )

        eloss, eloss_sigma = GetAnalyticalEloss(A, projectile_Z, Ein, material_name, length*0.1)
        dEcatima = rs.normal(loc=eloss, scale=eloss_sigma, size=10000)

        fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

        n_lines = 30
        ZMin = projectile_Z - 6
        ZMax = projectile_Z
        lines = []
        colors0 = []
        for i, h in enumerate(histories[0:n_lines]):
            for j in range(len(h) - 1):
                if h[j][1] > (total_length + length) * 0.1:
                    continue
                if h[j + 1][1] < total_length * 0.1:
                    continue
                lines.append([[h[j][1] * 10 - total_length, i], [h[j + 1][1] * 10 - total_length, i]])
                colors0.append(cm.viridis((h[j][0] - ZMin) / (ZMax - ZMin)))
        lc = mc.LineCollection(lines, colors=colors0, linewidths=3.8)

        ax = axes[0]
        ax.add_collection(lc)
        ax.autoscale()
        ax.set_xlim(0, length)
        ax.set_ylim(-int(n_lines * 0.05), n_lines)
        ax.set_yticks([0, 4, 9, 14, 19], [f"{i}" for i in [1, 5, 10, 15, 20]])
        ax.tick_params(labelleft=False, labelright=False)
        ax.set_xlabel("Thickness [mm]")
        ax.set_ylabel("Events")
        ax.set_title(material)

        ax = axes[1]
        for q in range(projectile_Z, projectile_Z - 4, -1):
            ax.plot(charges["length"] * 10, charges[q], label=f"Q={q}")
        ax.legend()
        ax.set_xlabel("Thickness [mm]")
        ax.set_xlim(0, None)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.4)
        ax.set_title("Charge state fraction")

        ax = axes[2]
        histRange = [np.mean(dEtotal) * 0.95, np.mean(dEtotal) * 1.05]
        for param, dE in zip(["dEcatima", "dEtotal (col+cc)"], [dEcatima, dEtotal]):
            label = f"{param}\nMean:{np.mean(dE):.3f}\nStdev:{np.std(dE):.3f} ({np.std(dE)/np.mean(dE):.2%})"
            ax.hist(dE, alpha=0.5, bins=50, range=histRange, label=label)
        ax.legend(loc="upper left", bbox_to_anchor=(0.8, 1))
        ax.set_xlabel("Energy loss [MeV/u]")
        ax.set_title("Energy loss distribution")
        ax.set_xlim(*histRange)
        ax.grid(alpha=0.4)
        st.pyplot(fig)

        st.write(f"Eout(simcc): {Ein - np.mean(dEtotal):.3f} MeV/u, Eloss(simcc): {np.mean(dEtotal):.3f} MeV/u, Eloss(CATIMA): {eloss:.3f} MeV/u ({(eloss)/(np.mean(dEtotal))-1:+.1%})")

        total_length += length
        Ein = Ein - np.mean(dEtotal)

    st.success("Calculation executed successfully!")
