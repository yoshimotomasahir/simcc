# SimCC
SimCCは物質中の電荷状態の変化をシミュレーションするためのPythonコードです。

## Dependencies

必須
- pycatima
- numpy

オプション
- jupyter
- streamlit
- matplotlib

## Features
- 電荷変化断面積や平均自由行程を計算
- MC法で物質中の電荷状態変化を計算
- 物質中のエネルギーロスを計算

## Getting Started
SimCCのインストールステップ

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

`error: externally-managed-environment` などと怒られたときは、 仮想環境をセットアップするか、`--break-system-packages`オプションをつけてインストールしてください。

## Usage

examples に Jupyter と streamlit 用のサンプルファイルあり。

streamlit を実行するには、

```sh
pip install streamlit
pip install matplotlib
streamlit run ./examples/streamlit_app.py
```
