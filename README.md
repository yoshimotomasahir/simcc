# SimCC
SimCC is a Python package for simulating charge-state changes in materials.

## Dependencies

Required:
- pycatima
- numpy

Optional:
- jupyter
- streamlit
- matplotlib

## Features
- Calculate charge-changing cross-sections and mean free paths
- Simulate charge-state changes in materials using the Monte Carlo method
- Calculate energy loss in materials

## Getting Started
Follow these steps to install SimCC:

1. Clone the repository:
    ```sh
    git clone https://ribfrepo.riken.jp/yoshimoto/simcc.git
    ```
2. Navigate to the project directory:
    ```sh
    cd simcc
    ```
3. Install using pip:
    ```sh
    pip install .
    ```
## Usage

Example files for Jupyter and Streamlit are available in the examples directory.

To run Streamlit, execute the following commands:

```sh
pip install streamlit matplotlib
streamlit run ./examples/app.py
```
