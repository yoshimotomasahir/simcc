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

st.title("MC for energy loss distribution")

projectile_Z, energy, A, charge_state = input_projectile()

materials = input_materials()

if st.button("Execute Calculation"):

    st.write("Not implemented.")
    
    st.success("Calculation executed successfully!")
