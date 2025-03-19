import streamlit as st
from uncertainties import ufloat
import plotly.graph_objects as go

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
    Energy = st.number_input("Energy", value=300, step=1)


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


def focus2int(FaFb):
    return int(FaFb.split("-")[0].replace("F", "")), int(FaFb.split("-")[1].replace("F", ""))


if "matrix" not in st.session_state:
    st.session_state.matrix = {}
for FaFb in ["F3-F5", "F7-F5"]:
    if FaFb not in st.session_state.matrix:
        st.session_state.matrix[FaFb] = get_matrix(*focus2int(FaFb))

for i, j in np.ndindex(X.shape):
    xval = X[i, j]
    F3X = ufloat(0, xval)
    F5X = ufloat(0, xval)
    F7X = ufloat(0, xval)
    F3A, F5A, F7A = ufloat(0, xval * 2), ufloat(0, xval * 2), ufloat(0, xval * 2)

    yval = Y[i, j]
    F3TOF = ufloat(0, yval / 1000)
    F7TOF = ufloat(0, yval / 1000)

    energy35 = Energy
    energy57 = Energy
    brho35 = energy2brho(energy35, A, Q)
    brho57 = energy2brho(energy57, A, Q)

    delta35, _ = calc_delta_ai(F3X, F5X, F5A, st.session_state.matrix["F3-F5"])
    delta57, _ = calc_delta_ai(F7X, F5X, F5A, st.session_state.matrix["F7-F5"])

    fl35 = calc_fl(F3X, F3A, delta35, st.session_state.matrix["F3-F5"], get_flc(3, 5))
    fl57 = calc_fl(F7X, F7A, delta57, st.session_state.matrix["F7-F5"], get_flc(5, 7))

    TOF35 = beta2tof(fl35, energy2beta(energy35))
    TOF35 += F3TOF
    TOF57 = beta2tof(fl57, energy2beta(energy57))
    TOF57 += F7TOF
    TOF37 = TOF35 + TOF57

    AOQ35, AOQ57, _, _, _, _ = calcAOQ_from_2brho(brho35 * (1 + delta35 * 0.01), brho57 * (1 + delta57 * 0.01), fl35, fl57, TOF37)
    AOQ35s[i, j] = AOQ35.s / AOQ35.n * 100
    AOQ57s[i, j] = AOQ57.s / AOQ35.n * 100

heatmap = go.Heatmap(z=AOQ35s, x=PPAC_error, y=TOF_error, colorscale="Viridis", opacity=0.7)
contour = go.Contour(z=AOQ35s, x=PPAC_error, y=TOF_error, colorscale="Blues", contours=dict(showlabels=True), contours_coloring="lines")
fig = go.Figure(data=[heatmap, contour])
fig.update_layout(title="BigRIPS A/Q resolution [%]", xaxis_title="PPAC resolution [mm]", yaxis_title="TOF resolution [ps]")
st.plotly_chart(fig)
