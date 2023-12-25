# pharmGKB-filter

Filter pharmGKB database by several rules

## dev

Prerequisites

- [Poetry](https://python-poetry.org/)
- [pyenv](https://github.com/pyenv/pyenv)

```bash
git clone https://github.com/Liyang-Zhang/pharmGKB-filter.git
cd pharmGKB-filter

pyenv install 3.11.6
poetry env use $(pyenv which python)

cat <<EOF > ".env"
export GPG_TTY=\$(tty)
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
source \$(poetry env info --path)/bin/activate
EOF

source .env
poetry install
pre-commit install
pre-commit run --all-files
```

## function

1. basic filter
2. add genome position
3. bed file filter

## usage

1. fill the path in the pharmgkb_filter/lung65_chemo.py file
2. ```bash
   python -m pharmgkb_filter.lung65_chemo
   ```
