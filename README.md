# SimCC
SimCCは、物質中の電荷状態の変化をシミュレーションするためのツールです

## Dependencies
- pycatima
- numpy

## Features
- 電荷変化断面積や平均自由行程を計算
- MC法で物質中の電荷状態変化を計算
- 物質中のエネルギーロスを計算

## Getting Started
To get started with SimCC, follow these steps:

1. Clone the repository:
    ```sh
    git clone https://ribfrepo.riken.jp/yoshimoto/simcc.git
    ```
2. Navigate to the project directory:
    ```sh
    cd simcc
    ```
3. pipでインストール:
    ```sh
    pip install .
    ```
4. Run the calculator:

## Usage

examples に Jupyter と streamlit 用のサンプルファイルあり。

```sh
pip install streamlit
pip install matplotlib
streamlit run ./examples/streamlit_app.py
```
