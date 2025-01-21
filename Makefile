.PHONY: venv install_maturin build develop python

BIN=venv/bin/
BIN_PIP=$(BIN)pip
BIN_MTN=$(BIN)maturin
PY_TOML=py/Cargo.toml

venv:
	python -m venv venv

install_maturin:
	$(BIN_PIP) install maturin

build:
	$(BIN_MTN) build -m $(PY_TOML) --release

develop:
	$(BIN_MTN) develop -m $(PY_TOML)

python:
	$(BIN)python
