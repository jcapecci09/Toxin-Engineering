## Conda enviornment

1. Ensure you are in the correct directory

```
cd Toxin-Engineering
```
2. Install conda enviornment

```
mamba env create -f mol-dyn.yml
```

3. Add environment as a jupyter kernel
```
python -m ipykernel install --user --name mol-dyn --display-name "Python (mol-dyn)"
```
