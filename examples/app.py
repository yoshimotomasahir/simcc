import streamlit as st

st.set_page_config(page_title="SimCC", page_icon="🌠")

st.title("Welcome to SimCC")
st.write("This web application allows you to simulate charge-states and energy loss of heavy ions in matter.")

st.write("## Features")

st.write("### Equilibrium state")
st.markdown("Displays the dependence of the charge states of heavy ions in equilibrium on the projectile atomic number, projectile energy, and target atomic number. [Go to Equilibrium state](./Equilibrium)")

st.write("### Transient state")
st.markdown("Calculates analytically charge state transitions of heavy ions in matter. Uses [pycatima](https://github.com/hrosiak/pycatima) to calculate the mean energy loss. Includes materials commonly used in the BigRIPS separator at RIBF. [Go to Transient state](./Transience)")

st.write("### Energy loss distribution")
st.markdown("Simulates energy loss event by event with Monte Carlo method. [Go to Energy loss distribution](./MC_Eloss_dist)")

st.write("## Repository Information")
st.markdown("The source code is available [here](https://github.com/yoshimotomasahir/simcc.git) managed by Masahiro Yoshimoto.")

st.sidebar.success("Select a page from the menu above.")
