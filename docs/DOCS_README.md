# Update Docs

```bash
cd docs
sphinx-apidoc -f -o source ../src/pyhiv
sphinx-autogen -o source/modules source/modules.rst
make clean
make html
```