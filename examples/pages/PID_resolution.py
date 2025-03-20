import streamlit as st
from uncertainties import ufloat
import uncertainties.umath as umath
import plotly.graph_objects as go
import plotly

import sys

sys.path.append("..")
from utils import *

st.set_page_config(page_title="PID resolution -SimCC-", page_icon="🌠")
st.header("PID resolution simulator for BigRIPS")
st.write("Under development.")

col1, col2, col3, col4 = st.columns(4)
with col1:
    Z = st.number_input("Atomic Number (Z)", value=50, step=1)
with col2:
    A = st.number_input("Mass Number (A)", value=132, step=1)
with col3:
    Q = st.number_input("Charge State (Q)", value=50, step=1)
with col4:
    Energy = st.number_input("Energy [MeV/u]", value=300, step=1)

col1, col2, col3 = st.columns(3)
with col1:
    deg_list = [
        "Al 0.7 mm 0.61 mrad",
        "Al 1.0 mm 0.43 mrad",
        "Al 1.0 mm 0.95 mrad",
        "Al 1.5 mm 1.187 mrad",
        "Al 2.0 mm 1.60 mrad",
        "Al 2.2 mm 1.80 mrad",
        "Al 2.5 mm 2.14 mrad",
        "Al 3.0 mm 2.62 mrad",
        "Al 3.5 mm 2.80 mrad",
        "Al 3.5 mm 3.00 mrad",
        "Al 4.5 mm 3.70 mrad",
        "Al 5.0 mm 4.25 mrad",
        "Al 7.0 mm 5.969 mrad",
        "Al 8.0 mm 7.3129 mrad",
        "Al 9.0 mm 7.714 mrad",
        "Al 10.0 mm 8.587 mrad",
    ]
    F5_deg = st.selectbox("F5 degrader", options=deg_list, index=4)
    F5_deg_center = float(F5_deg.split()[1]) * 0.1  # cm
    F5_deg_angle = float(F5_deg.split()[3])  # mrad
with col2:
    F5_deg_error = st.selectbox(
        "Thickness ununiformity",
        options=[
            "0 μm",
            "1 μm",
            "2 μm",
            "5 μm",
            "10 μm",
            "15 μm",
            "20 μm",
            "25 μm",
        ],
        index=4,
    )
    F5_deg_error = float(F5_deg_error.split()[0]) * 0.1 * 0.001
with col3:
    Straggling_enhancement = st.number_input("Straggling enhancement by cc", value=1.0, step=0.1, min_value=1.0)


eloss, straggling = GetAnalyticalEloss(A, Z, Energy, "Al", F5_deg_center)
F5_deg_thickness = ufloat(F5_deg_center, F5_deg_error)
F5_deg_thickness *= ufloat(1, straggling / eloss * Straggling_enhancement)

energy35_nominal = Energy
energy57_nominal = Energy - eloss
brho35_nominal = energy2brho(energy35_nominal, A, Q)
brho57_nominal = energy2brho(energy57_nominal, A, Q)


def calc_delta_ai(Xi, Xf, Af, mat):
    # mat Fi -> Ff
    D = mat[1][1] * mat[5][0] - mat[1][0] * mat[5][1]
    delta = mat[1][1] * Xf - mat[1][0] * Af - (mat[1][1] * mat[0][0] - mat[1][0] * mat[0][1]) * Xi
    delta = delta / D
    recoAi = mat[5][1] * Xf - mat[5][0] * Af - (mat[5][1] * mat[0][0] - mat[5][0] * mat[0][1]) * Xi
    recoAi = -recoAi / D
    return delta, recoAi


def get_matrix(Fa, Fb):
    import requests

    url = "https://ribf.riken.jp/BigRIPSInfo/optics/fig/bStandard.txt"
    response = requests.get(url)
    assert response.status_code == 200
    data = response.text
    lines = data.splitlines()
    mat = []
    for i, line in enumerate(lines):
        if line.replace("\n", "").replace("\r", "") == "F{}-F{}".format(Fa, Fb):
            for j in range(6):
                mat.append([float(f) for f in lines[i + j + 1].split()[0:5]])
            assert lines[i + 7][0] == "-"
            return mat


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


def get_flc(Fa, Fb):
    if Fa == 3 and Fb == 5:
        return 23488.0
    elif Fa == 5 and Fb == 7:
        return 23488.0
    elif Fa == 3 and Fb == 7:
        return 23488.0 * 2
    else:
        raise


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


PPAC_error = np.linspace(0.1, 1.0, 19)
TOF_error = np.linspace(10, 120, 12)
X, Y = np.meshgrid(PPAC_error, TOF_error)
AOQ35s = np.zeros_like(X)
AOQ57s = np.zeros_like(X)
Zdegs = np.zeros_like(X)


def focus2int(FaFb):
    return int(FaFb.split("-")[0].replace("F", "")), int(FaFb.split("-")[1].replace("F", ""))


if "matrix" not in st.session_state:
    st.session_state.matrix = {}
for FaFb in ["F3-F5", "F7-F5"]:
    if FaFb not in st.session_state.matrix:
        st.session_state.matrix[FaFb] = get_matrix(*focus2int(FaFb))


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


for i, j in np.ndindex(X.shape):
    xval = X[i, j]
    F3X = ufloat(0, xval)
    F5X = ufloat(0, xval)
    F7X = ufloat(0, xval)
    F3A, F5A, F7A = ufloat(0, xval * 2), ufloat(0, xval * 2), ufloat(0, xval * 2)

    yval = Y[i, j]
    F3TOF = ufloat(0, yval / 1000)
    F7TOF = ufloat(0, yval / 1000)

    delta35, _ = calc_delta_ai(F3X, F5X, F5A, st.session_state.matrix["F3-F5"])
    delta57, _ = calc_delta_ai(F7X, F5X, F5A, st.session_state.matrix["F7-F5"])

    fl35 = calc_fl(F3X, F3A, delta35, st.session_state.matrix["F3-F5"], get_flc(3, 5))
    fl57 = calc_fl(F7X, F7A, delta57, st.session_state.matrix["F7-F5"], get_flc(5, 7))

    TOF35 = beta2tof(fl35, energy2beta(energy35_nominal))
    TOF35 += F3TOF
    TOF57 = beta2tof(fl57, energy2beta(energy57_nominal))
    TOF57 += F7TOF
    TOF37 = TOF35 + TOF57

    brho35 = brho35_nominal * (1 + delta35 * 0.01)
    brho57 = brho57_nominal * (1 + delta57 * 0.01)
    AOQ35, AOQ57, beta35, beta57, _, _ = calcAOQ_from_2brho(brho35, brho57, fl35, fl57, TOF37)
    Zdeg = calc_zdeg(F5X, F5A, beta35, beta57, brho35, brho57, F5_deg_thickness, F5_deg_angle, Z)

    AOQ35s[i, j] = AOQ35.s / AOQ35.n * 100
    AOQ57s[i, j] = AOQ57.s / AOQ35.n * 100
    Zdegs[i, j] = Zdeg.s

heatmap1 = go.Heatmap(z=AOQ35s, x=PPAC_error, y=TOF_error, colorscale="Viridis", opacity=0.7)
contour1 = go.Contour(z=AOQ35s, x=PPAC_error, y=TOF_error, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines")
fig = go.Figure(data=[heatmap1, contour1])
fig.update_layout(title="A/Q resolution [%] (std.dev.)", xaxis_title="PPAC resolution [mm]", yaxis_title="TOF resolution [ps]",margin=dict(l=5, r=5, t=30, b=5),width=1000,height=300)
st.plotly_chart(fig)

heatmap2 = go.Heatmap(z=Zdegs, x=PPAC_error, y=TOF_error, colorscale="Viridis", opacity=0.7)
contour2 = go.Contour(z=Zdegs, x=PPAC_error, y=TOF_error, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines")
fig = go.Figure(data=[heatmap2, contour2])
fig.update_layout(title="Zdeg resolution (std.dev.)", xaxis_title="PPAC resolution [mm]", yaxis_title="TOF resolution [ps]",margin=dict(l=5, r=5, t=30, b=5),width=1000,height=300)
st.plotly_chart(fig)
