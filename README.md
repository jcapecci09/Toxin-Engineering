#  Toxin-Engineering

This project aims to engineer the scorpion toxin LQHIII to bind cardiac proteins. The workflow uses the experimentally verified NaV1.5–LQHIII complex as a starting point, modifies NaV1.5 to mimic cardiac proteins, and then engineers LQHIII to bind the resulting cardiac-like complex.

## ✨ Highlights

This project provides two reusable PyRosetta-based tools for protein engineering and structural analysis:

- `insert_mutation`: Replaces one sequence motif with another motif in a protein structure.
- `AffinityOptimizer`: Redesigns a protein–protein binding interface to identify mutations that may improve binding affinity.

**Note:** The underlying implementation of these tools can be found in `utils.py`, while examples of their usage can be found in `demo.ipynb`.

The complete analysis is detailed in `project.ipynb` and consists of the following components:

1. Explore the experimentally verified wildtype NaV1.5–LQHIII complex.
2. Use `insert_mutation` to redesign NaV1.5 to mimic cardiac proteins.
3. Use `AffinityOptimizer` to engineer LQHIII and generate variants with potentially improved binding affinity for the cardiac-like NaV1.5.
4. Explore and evaluate the resulting variants.

## 🧬 Overview

   
## 🛠️ Methods

This section provides a brief overview of the methods used in this analysis. The workflow consisted of two major phases, described below:

### Phase 1: Generating a Cardiac-Like NaV1.5

The goal of Phase 1 was to mimic cardiac protein sequences by sequentially introducing mutations into NaV1.5. Mutations were introduced one amino acid at a time until the **ABLIM1 MSSSP motif**, associated with GSK3β recognition, was incorporated.

**NaV1.5:** `GTVLSDIIQKY` → **Cardiac-like:** `AQPMSSSPKET`

```python
mutant_pose = insert_mutation('pdb_files/7K18_relaxed.pdb', 'GTVLSDIIQKY', 'AQPMSSSPKET')
```

### Phase 2: Generate LQHIII Variants

The goal of Phase 2 was to mutate **LQHIII** to improve its predicted binding to the cardiac-like NaV1.5. A **simulated annealing** approach was used to explore mutations and optimize binding scores.

- Randomly introduce mutations and relax the surrounding structure
- Score each variant using **distance-based interactions** and **ΔG**
- Accept favorable mutations and occasionally accept unfavorable ones to explore alternative sequences
- Gradually lower the temperature to favor better-performing variants
- Repeat for a set number of iterations

**Scoring:**
- **Distance Score:** Evaluates targeted interactions such as hydrogen bonds and salt bridges.
- **ΔG:** Evaluates overall predicted binding energetics, with more negative values indicating stronger predicted binding.

The baseline hyperparameters are shown below:

```python
# set amino acids to allow mutation
aa_mutate = [107, 108, 109, 110, 113, 114, 117, 118, 120, 121,
                141, 142, 157, 159, 160, 161, 162, 166]
opt1 = AffinityOptimizer(
    'pdb_files/7K18_AQPMSSSPKET_mutant.pdb',
    scoring_function='DDG',
    pos_to_mutate=aa_mutate,
    temp=20,
    cooling_rate=0.95,
    relax=True,
    relax_every=1,
    early_stop_iter=75,
    number_steps=500,
    quiet=True
)
```

**Note:** A total of 160 variants were generated: 80 using the **Distance Score** and 80 using the **ΔG Score**. Hyperparameters such as `temp`, `cooling_rate`, and `relax` were adjusted between runs to explore different optimization conditions.

## 📊 Results 

### Runtime
![](Data/figures/runtime.png)

### Best Variants
![](Data/figures/top_variants_dg_dis_scatter_wt.png)
![](Data/figures/Best_variants_wt_scatter.png)
![](Data/figures/top_variants_sequences.png)

## ⚙️ Usage Instructions

### <u>Project Setup</u>

1. Clone the repository and navigate into the project directory:
   
```
git clone https://github.com/jcapecci09/Toxin-Engineering.git
cd Toxin-Engineering
```

### <u>Environment Setup</u>

This analysis requires a Mamba environment. Ensure you have [Mamba installed](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) before proceeding.

1. Create the environment from the file:

```
mamba env create -f mol-dyn.yml
```

2. Activate the environment:
   
```
mamba activate mol-dyn
```

3. Register the environment as a Jupyter kernel:

```
python -m ipykernel install --user --name mol-dyn --display-name "Python (mol-dyn)"
```

4. Launch Jupyter Notebook or VS Code and select **Python (mol-dyn)** as your active kernel.

5. Work through `project.ipynb` to produce results


## 👤 Author

I'm [Jimmy Capecci](https://github.com/jcapecci09), a Bioinformatics graduate student at Loyola University Chicago. I completed this project during my internship at the Stritch School of Medicine in the Peter Kekenes-Huskey Lab, where I explored computational approaches to protein engineering and structural analysis.

## 🙏 Acknowledgements

I'd like to thank Dr. Kekenes-Huskey and Alec Loftus at the Stritch School of Medicine for their guidance and support throughout this project.
