#!/bin/bash

python3 -m venv scanpy_env
source scanpy_env/bin/activate

pip install --upgrade pip

pip install \
    scanpy \
    matplotlib \
    seaborn \
    numpy \
    scipy \
    pandas

echo "Installation complete."
echo "Activate with: source scanpy_env/bin/activate"
