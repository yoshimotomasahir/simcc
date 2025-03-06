import streamlit as st

st.set_page_config(page_title="SimCC", page_icon="🌠")

st.title("Welcome to SimCC")
st.write("This application allows you to simulate charge-state changes in materials.")

st.write("## Repository Information")
repo_url = "https://github.com/yoshimotomasahir/simcc.git"
message = f"The source code is available [here]({repo_url})."
st.write(message)

st.sidebar.success("Select a page from the menu above.")
