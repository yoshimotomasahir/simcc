from setuptools import setup, find_packages

setup(
    name="simcc",
    version="0.1.0",
    author="Masahiro Yoshimoto",
    author_email="masahiro.yoshimoto@riken.jp",
    description="A package for simulating charge-state changing.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yoshimotomasahir/simcc",
    packages=find_packages(),
    package_data={"simcc": ["ChargeStates/*.txt"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=["numpy>=1.26.0", "pycatima>=1.96"],
)
