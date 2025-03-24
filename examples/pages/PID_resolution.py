import streamlit as st
from simcc import *
from uncertainties import ufloat
import uncertainties.umath as umath
import plotly.graph_objects as go
import plotly

import sys

sys.path.append("..")
from utils import *

st.set_page_config(page_title="PID Resolution -SimCC-", page_icon="🌠")
st.header("PID Resolution Simulator for BigRIPS")
st.write("Energy loss other than the F5 degrader is neglected. The resolution is the **standard deviation**.")
st.markdown("The source code is available on [Github](https://github.com/yoshimotomasahir/simcc/blob/main/examples/pages/PID_resolution.py).")

col1, col2, col3, col4 = st.columns(4)
with col1:
    Z = st.number_input("Atomic Number (Z)", value=50, step=1, min_value=1)
with col2:
    A = st.number_input("Mass Number (A)", value=132, step=1, min_value=Z)
with col3:
    dQ = st.number_input("Z - Charge State (Q)", value=0, step=1, max_value=4, min_value=0)
    Q = Z - dQ
with col4:
    Energy = st.number_input("Energy [MeV/u] at F3", value=300, step=1, min_value=50)

col1, col2, col3 = st.columns(3)
with col1:
    F5_deg = st.selectbox("F5 Degrader", options=F5_deg_list, index=5)
    F5_deg_center = float(F5_deg.split()[1]) * 0.1  # cm
    F5_deg_angle = float(F5_deg.split()[3])  # mrad
with col2:
    F5_deg_error = st.number_input("Thickness Uniformity [μm]", value=10, step=1, min_value=0)  # um
    F5_deg_error = F5_deg_error * 0.1 * 0.001  # cm
with col3:
    Straggling_enhancement = st.number_input("Eloss Straggling Enhancement by cc", value=1.0, step=0.1, min_value=1.0)


def calc_delta_ai(Xi, Xf, Af, mat):
    # mat Fi -> Ff
    D = mat[1][1] * mat[5][0] - mat[1][0] * mat[5][1]
    delta = mat[1][1] * Xf - mat[1][0] * Af - (mat[1][1] * mat[0][0] - mat[1][0] * mat[0][1]) * Xi
    delta = delta / D
    recoAi = mat[5][1] * Xf - mat[5][0] * Af - (mat[5][1] * mat[0][0] - mat[5][0] * mat[0][1]) * Xi
    recoAi = -recoAi / D
    return delta, recoAi


def calc_fl(Xi, Ai, delta, mat, flc):
    # mat Fi -> Ff
    fl = flc
    assert flc > 0
    assert len(mat) == 6
    fl += (mat[0][4] * Xi + mat[1][4] * Ai + mat[5][4] * delta) * (-1.0)
    return fl


def calcAOQ(Brho, beta):
    gamma = 1.0 / ((1 - beta) * (1 + beta)) ** 0.5
    AOQ = Brho * clight / mnucleon / beta / gamma
    return AOQ


def calcAOQ_from_2brho(brho35, brho57, fl35, fl57, TOF37):
    a = brho57 / brho35
    b = (a**2 * clight**2 * TOF37**2 + a**2 * (a + 1) * (a - 1) * fl35**2 + (1 - a) * (1 + a) * fl57**2) ** 0.5

    beta35 = b * fl35 + fl57 * clight * TOF37
    beta35 = beta35 / (b * clight * TOF37 + (1 - a) * (1 + a) * fl35 * fl57)

    beta57 = b * fl35 + fl57 * clight * TOF37
    beta57 = beta57 / (clight**2 * TOF37**2 + (a + 1) * (a - 1) * fl35**2)

    gamma35 = 1.0 / ((1 - beta35) * (1 + beta35)) ** 0.5
    gamma57 = 1.0 / ((1 - beta57) * (1 + beta57)) ** 0.5
    AOQ35 = brho35 * clight / mnucleon / beta35 / gamma35
    AOQ57 = brho57 * clight / mnucleon / beta57 / gamma57
    return AOQ35, AOQ57, beta35, beta57, gamma35, gamma57


def calc_zdeg(F5X, F5A, beta35, beta57, brho35, brho57, d0, angle, Z):
    zpos = 0
    pos = F5X + zpos * umath.tan(F5A / 1000)
    t = d0 + pos * angle / 1000
    t = t / umath.cos(F5A / 1000)
    ionpair = 6156  # Al
    de_v = t * (umath.log(ionpair * beta35 * beta35) - umath.log(1 - beta35 * beta35) - beta35 * beta35)
    delE = brho35 / beta35 - brho57 / beta57
    a = Z / ((delE / de_v) * beta35 * beta35)
    zdeg = a.n * (delE / de_v) * beta35 * beta35
    return zdeg


if F5_deg_center > 0:
    eloss, straggling = GetAnalyticalEloss(A, Z, Energy, "Al", F5_deg_center)
    F5_deg_thickness = ufloat(F5_deg_center, F5_deg_error)
    F5_deg_thickness *= ufloat(1, straggling / eloss * Straggling_enhancement)
else:
    eloss = 0
    F5_deg_thickness = ufloat(0, 0)


energy35_nominal = Energy
energy57_nominal = Energy - eloss
brho35_nominal = energy2brho(energy35_nominal, A, Q)
brho57_nominal = energy2brho(energy57_nominal, A, Q)
if "matrixF3F5" not in st.session_state:
    st.session_state.matrixF3F5 = get_matrix(3, 5)
if "matrixF7F5" not in st.session_state:
    st.session_state.matrixF7F5 = get_matrix(7, 5)
matrix35 = st.session_state.matrixF3F5
matrix75 = st.session_state.matrixF7F5

position_errors = np.linspace(0.0, 0.8, 17)
timing_errors = np.linspace(0, 90, 10)
mesh_grid = np.meshgrid(position_errors, timing_errors)
relative_AOQ35s = np.zeros_like(mesh_grid[0])
absolute_AOQ35s = np.zeros_like(mesh_grid[0])
Zdegs = np.zeros_like(mesh_grid[0])


for i, j in np.ndindex(mesh_grid[0].shape):
    position_error = mesh_grid[0][i, j]
    F3X, F5X, F7X = ufloat(0, position_error), ufloat(0, position_error), ufloat(0, position_error)
    F3A, F5A, F7A = ufloat(0, position_error * 2), ufloat(0, position_error * 2), ufloat(0, position_error * 2)

    timing_error = mesh_grid[1][i, j]
    timingF3 = ufloat(0, timing_error / 1000)
    timingF7 = ufloat(0, timing_error / 1000)

    delta35, _ = calc_delta_ai(F3X, F5X, F5A, matrix35)
    delta57, _ = calc_delta_ai(F7X, F5X, F5A, matrix75)

    fl35 = calc_fl(F3X, F3A, delta35, matrix35, 23488.0)
    fl57 = calc_fl(F7X, F7A, delta57, matrix75, 23488.0)

    TOF35 = beta2tof(fl35, energy2beta(energy35_nominal))
    TOF35 += timingF3
    TOF57 = beta2tof(fl57, energy2beta(energy57_nominal))
    TOF57 += timingF7
    TOF37 = TOF35 + TOF57

    brho35 = brho35_nominal * (1 + delta35 * 0.01)
    brho57 = brho57_nominal * (1 + delta57 * 0.01)
    AOQ35, AOQ57, beta35, beta57, _, _ = calcAOQ_from_2brho(brho35, brho57, fl35, fl57, TOF37)
    relative_AOQ35s[i, j] = AOQ35.s / AOQ35.n * 100
    absolute_AOQ35s[i, j] = AOQ35.s

    if F5_deg_thickness.n > 0:
        Zdeg = calc_zdeg(F5X, F5A, beta35, beta57, brho35, brho57, F5_deg_thickness, F5_deg_angle, Z)
        Zdegs[i, j] = Zdeg.s


heatmap = go.Heatmap(z=relative_AOQ35s, x=position_errors, y=timing_errors, colorscale="Viridis", opacity=0.7, name="")
contour = go.Contour(z=relative_AOQ35s, x=position_errors, y=timing_errors, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines", name="")
fig = go.Figure(data=[heatmap, contour])
fig.update_layout(title="Relative A/Q Resolution [%]", xaxis_title="PPAC Position resolution [mm]", yaxis_title="Timing resolution [ps]", margin=dict(l=5, r=5, t=30, b=5), width=600, height=300)
st.plotly_chart(fig)

if F5_deg_thickness.n > 0:
    heatmap = go.Heatmap(z=Zdegs, x=position_errors, y=timing_errors, colorscale="Viridis", opacity=0.7, name="")
    contour = go.Contour(z=Zdegs, x=position_errors, y=timing_errors, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines", name="")
    fig = go.Figure(data=[heatmap, contour])
    fig.update_layout(title="Absolute Zdeg Resolution", xaxis_title="PPAC Position resolution [mm]", yaxis_title="Timing resolution [ps]", margin=dict(l=5, r=5, t=30, b=5), width=600, height=300)
    st.plotly_chart(fig)

heatmap = go.Heatmap(z=absolute_AOQ35s * 1000, x=position_errors, y=timing_errors, colorscale="Viridis", opacity=0.7, name="")
contour = go.Contour(z=absolute_AOQ35s * 1000, x=position_errors, y=timing_errors, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines", name="")
fig = go.Figure(data=[heatmap, contour])
fig.update_layout(title="Absolute A/Q Resolution [10^-3]", xaxis_title="PPAC Position resolution [mm]", yaxis_title="Timing resolution [ps]", margin=dict(l=5, r=5, t=30, b=5), width=600, height=300)
st.plotly_chart(fig)
