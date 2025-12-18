import streamlit as st
import numpy as np
import pycatima as catima
from simcc import *
import sys

sys.path.append("..")
from utils import *

import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd

st.set_page_config(page_title="Material Optimisation -SimCC-", page_icon="🌠")
st.header("Material Optimisation")

projectile_Z, energy, A, charge_state = input_projectile()

dE = st.number_input(
    "Energy loss dE [MeV/u]",
    value=1.0,
    min_value=0.1,
    format="%.3f",
    step=1.0,
)


def GetThickness(A, Z, energy, dE, material, tmax=200.0):
    if energy - dE >= energy:
        return 0.0
    if energy - dE <= 0:
        return np.nan

    result = GetMaterial(material)
    zts, m_fractions = result["zts"], result["m_fractions"]

    def GetCAtimaCompound(zts, m_fractions):
        return [[0, z, m] for m, z in zip(m_fractions, zts)]

    compound = GetCAtimaCompound(zts, m_fractions)
    mat = catima.Material(compound)
    config = catima.Config()
    config.z_effective = 1

    def f(t):
        mat.thickness(t)  # t: g/cm^2
        res = catima.calculate(catima.Projectile(A=A, Z=Z, Q=0, T=energy), mat, config)
        return res.Eout - (energy - dE)  # 0 になればOK

    # bracket を作る（ここだけ最小限）
    t_hi = 1e-4
    while t_hi < tmax and f(t_hi) > 0:
        t_hi *= 1.6
    if t_hi >= tmax:
        return np.nan

    from scipy.optimize import root_scalar

    sol = root_scalar(f, bracket=(0.0, t_hi), method="brentq")
    return sol.root if sol.converged else np.nan


Probs = []
Thicknesses = []

materials = list(range(1, 93))
for material in materials:
    thickness = GetThickness(A, projectile_Z, energy, dE, material)
    Thicknesses.append(thickness * 1000)

    MFP = GetMFP(zp=projectile_Z, energy=energy + dE / 2, material=material)
    P0 = np.zeros(7)
    P0[charge_state] = 1
    P0 = GetAnalyticalProb(MFP, thickness, charge_state=P0)
    Probs.append(P0)

fig = make_subplots(subplot_titles=["Charge-state probability"])
for j in range(7):
    fig.add_trace(go.Scatter(x=materials, y=np.array(Probs).T[j].T, mode="lines", name=f"Z-Q={j}"))
fig.update_layout(xaxis_title="Stripper Z")
st.plotly_chart(fig)

Elements = ["Be", "Carbon", "Al", "Cu", "Nb", "Ta"]
rows = []
for element in Elements:
    try:
        result = GetMaterial(element)
        zts, density = result["zts"], result["density"]
        k = zts[0] - 1
        thickness = Thicknesses[k]
        rows.append(
            {
                "Stripper": element,
                "Z": zts[0],
                "Thickness [um]": thickness / density / 100 * 1000,
                "Full-strip [%]": Probs[k][0] * 100,
                "H-like [%]": Probs[k][1] * 100,
                "He-like [%]": Probs[k][2] * 100,
                "Li-like [%]": Probs[k][3] * 100,
                "x [mg/cm2]": thickness,
                "ρ [mg/cm3]": density * 1000,
            }
        )
    except:
        pass

df = pd.DataFrame(rows).set_index("Stripper")
st.write("Equivalent stripper thickness and charge states for a given energy loss.")
st.dataframe(df.round(1), width="stretch")

fig = make_subplots(subplot_titles=["Equivalent thickness(x) [mg/cm2]"])
fig.add_trace(go.Scatter(x=materials, y=np.array(Thicknesses), mode="lines"))
fig.update_layout(xaxis_title="Stripper Z")
st.plotly_chart(fig)

st.success("Calculation executed successfully!")
