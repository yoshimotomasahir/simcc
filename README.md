# SimCC
SimCC is a Python package for simulating charge-state changes in materials.

You can try this package at https://simcc-web.streamlit.app/ .

## Dependencies

Required:
- pycatima
- numpy
- scipy

Optional:
- jupyter
- streamlit
- matplotlib
- plotly

## Features
- Calculate charge-changing cross-sections and mean free paths
- Simulate charge-state changes in materials using the Monte Carlo method
- Calculate energy loss in materials

## Getting Started
Follow these steps to install SimCC:

1. Create and activate a virtual environment:
    ```sh
    sudo apt update
    sudo apt install python3-venv python3-dev
    cd ~/venv
    python3 -m venv venv
    source venv/bin/activate  # On macOS/Linux
    venv\Scripts\activate  # On Windows
    ```
2. Install `build-essential` and `cmake`:
    ```sh
    sudo apt install build-essential cmake
    ```
3. Clone the repository:
    ```sh
    cd ~
    git clone https://github.com/yoshimotomasahir/simcc.git
    ```
4. Navigate to the project directory:
    ```sh
    cd simcc
    ```
5. Install using pip:
    ```sh
    pip install .
    ```
## Usage

Example files for Jupyter and Streamlit are available in the examples directory.

To run Streamlit, execute the following commands:

```sh
streamlit run ./examples/app.py
```

## Notes

Note that this simulation does not consider transitions of two-electrons capture and two-electrons removal.

Note that the cross sections does not match perfectly with the [GLOBAL](https://web-docs.gsi.de/~weick/charge_states/) ones, see the [cross section data](/simcc/ChargeStates).
