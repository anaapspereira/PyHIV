# How to update documentation:

## 1. From the root of the documentation folder (docs/), run:

```bash
sphinx-apidoc -f -e -M -o source ../src/pyhiv
```

## 2. Then, build the documentation:

```bash
make clean
make html
```
