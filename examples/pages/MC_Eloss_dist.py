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
st.header("MC for energy loss distribution")

projectile_Z, energy, A, charge_state = input_projectile()

materials = input_materials()

if st.button("Execute Calculation"):

    expanded_materials = get_expanded_materials(materials)
    Ein = energy
    histories = None
    total_length = 0
    rs = np.random.RandomState(1)
    rs_mc = np.random.RandomState(1)
    for k, material in enumerate(expanded_materials):
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
        elossm1, _ = GetAnalyticalEloss(A, projectile_Z-1, Ein, material_name, length*0.1)
        N=10000
        dEcatima = rs.normal(loc=eloss, scale=eloss_sigma, size=N)

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

        st.markdown("""
|Ein [MeV/u]| |Eout|Eloss|Sigma|Zres<sup>*1</sup>|Eres<sup>*2</sup>|
|--|--|--|--|--|--|--|
|{}|{}|{}|{}|{}|{}|{}|
|{}|{}|{}|{}|{}|{}|{}|
|{}|{}|{}|{}|{}|{}|{}|
""".format(
f"{Ein:.6g}",
"CATIMA",
f"{Ein - eloss:.3f}",
f"{eloss:.3f}",
f"{eloss_sigma:.4f}",
f"{(eloss-elossm1)/eloss_sigma:.3f}σ",
f"{eloss_sigma/eloss:.3%}",

f"<sup>{A}</sup>{z2symbol[projectile_Z]}"+(f"<sup>{projectile_Z - charge_state}+</sup>" if k==0 else "<sup>all+</sup>"),
"simcc",
f"{Ein - np.mean(dEtotal):.3f}",
f"{np.mean(dEtotal):.3f}",
f"{np.std(dEtotal):.4f}",
f"{(eloss-elossm1)/np.std(dEtotal):.3f}σ",
f"{np.std(dEtotal)/np.mean(dEtotal):.3%}",

f"{material_name} {length:.5g} mm",
"cc/ATIMA",
f"100{(Ein - np.mean(dEtotal))/(Ein - eloss)-1:+.4%}",
f"100{(np.mean(dEtotal))/(eloss)-1:+.3%}",
f"{(np.std(dEtotal))/(eloss_sigma):.2f}<sup>*3</sup>",
f"{(eloss-elossm1):.4f}<sup>*4</sup>",
f"{(eloss-elossm1)/eloss:.3%}<sup>*5</sup>",
),unsafe_allow_html=True)

        csv_text="dEcatima\tdEtotal\n"+"\n".join(f"{a}\t{b}" for a,b in zip(dEcatima,dEtotal))
        with st.expander("Show raw data as TSV format",expanded=False):
            st.code(csv_text,language="text")
        
        total_length += length
        Ein = Ein - np.mean(dEtotal)

    st.success(f"Calculation (N={N}) executed successfully!")

st.html("""
*1: Sigma / (Eloss<sub>Z</sub> - Eloss<sub>Z-1</sub>)<br>
*2: Sigma / Eloss<br>
*3: Eloss Straggling Enhancement by cc<br>
*4: Eloss<sub>Z</sub> - Eloss<sub>Z-1</sub><br>
*5: (Eloss<sub>Z</sub> - Eloss<sub>Z-1</sub>)/Eloss<sub>Z</sub><br>
""")
