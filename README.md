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

Intrinsically disordered regions (IDRs) are flexible regions of proteins that do not adopt a single stable structure, making them difficult to study using traditional structural methods (Sun et al., 2025). Despite this challenge, IDRs play important roles in regulating protein interactions and often contain **phosphorylation sites** that can alter protein function.

**ABLIM1**, a cardiac myofilament protein, contains IDRs with phosphorylation sites that may regulate its interactions with other cardiac proteins. This project uses the **ABLIM1 MSSSP motif**, a phosphorylation-associated sequence, as a basis for engineering a cardiac-like **NaV1.5** binding region.

**LQHIII**, a disulfide-rich scorpion toxin that targets voltage-gated sodium channels, was selected as the protein engineering target due to its high structural stability and functional activity under physiological conditions (He et al, 2025). The engineered NaV1.5 structure was then used to optimize **LQHIII** variants with **PyRosetta** and simulated annealing to explore changes in predicted protein–protein binding.

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

**Note**: Native LQHIII–NaV1.5 interactions were examined to identify key contacts near the engineered binding region. For the cardiac-like complex, K64–E1616, H15–S1612, and H43–S1611 were selected as potential salt bridge and hydrogen bond targets for distance-based optimization. The full interaction analysis is provided in project.ipynb.

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
**Figure 1. Runtime of `AffinityOptimizer`.** Runtime depends on the size of the complex, selected hyperparameters, and computational system. Using the baseline hyperparameters described above and the truncated 167-AA complex, `AffinityOptimizer` typically required **5–10 minutes per variant**. The **Distance Score** generally ran faster than the **ΔG Score**. Disabling relaxation significantly reduced runtime but may reduce structural quality, as the model is less constrained to maintain a realistic complex.

### Best Variants
![](Data/figures/top_variants_dg_dis_scatter_wt.png)
**Figure 2. Top 15 variants from Distance Score and top 15 variants from ΔG Score compared with the Wildtype and Initial structures.** The **Wildtype** represents the experimentally validated complex before Phase 1, while the **Initial** structure represents the cardiac-like NaV1.5 complex generated after Phase 1. Variants optimized using the **Distance Score** generally showed improvements in both targeted interaction distances and overall predicted binding energy (ΔG). In contrast, variants optimized directly for **ΔG** achieved more favorable binding energies but generally had higher Distance Squared Error (DSE), indicating that improved overall binding energetics did not necessarily correspond to improved targeted interactions.  

<br>

![](Data/figures/Best_variants_wt_scatter.png)
**Figure 3. Top candidate variants based on Distance Squared Error (DSE) and Rosetta estimated ΔG.** Variants with DSE values below 7 and ΔG values below −24.5 were selected as potential candidates for binding to the cardiac-like NaV1.5. The **Wildtype** is shown as a reference for comparison. Several hyperparameter conditions produced variants with scores similar to the Wildtype. The **Distance Score optimization with relaxation disabled** produced the lowest DSE while maintaining a comparable ΔG, suggesting that targeted interaction distances could be optimized without substantially reducing predicted binding affinity.  

<br>

![](Data/figures/top_variants_sequences.png)
**Figure 4. Sequence Comparison of Top Candidate Variants.** Amino acid sequences of the top candidate variants compared with the relaxed wildtype reference. Gray indicates sequence differences, while red indicates residues involved in important interactions.  

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

## 🧾 Conclusion
This project explored how computational protein engineering can be used to redesign the LQHIII scorpion toxin toward a cardiac-like NaV1.5 binding environment. By first modifying NaV1.5 to mimic the ABLIM1 sequence and then using AffinityOptimizer to explore mutations in LQHIII, I was able to generate and compare variants across different optimization conditions.

One of the main takeaways from the results was that optimizing for targeted interactions and optimizing for overall binding energy did not always produce the same variants. The Distance Score produced candidates with better interaction geometry and, in several cases, maintained ΔG values comparable to the reference structures. This suggests that considering both structural interactions and overall binding energetics can be useful when selecting candidates.

The variants generated here are computational predictions and still need to be experimentally tested to determine whether they actually improve binding to cardiac proteins. However, this project provides a starting point for using PyRosetta and simulated annealing to explore protein–protein interface engineering and prioritize candidates for future experimental work.

## 👤 Author

I'm [Jimmy Capecci](https://github.com/jcapecci09), a Bioinformatics graduate student at Loyola University Chicago. I completed this project during my internship at the Stritch School of Medicine in the Peter Kekenes-Huskey Lab, where I explored computational approaches to protein engineering and structural analysis.

## 📚 Sources

He, D., Lei, Y., Qin, H., Cao, Z., & Kwok, H. F. (2025). Deciphering scorpion toxin-induced pain: Molecular mechanisms and ion channel dynamics. *International Journal of Biological Sciences, 21*(7), 2921–2934. [https://doi.org/10.7150/ijbs.109713](https://doi.org/10.7150/ijbs.109713)

Sun, B., Loftus, A., Beh Goh Beh, B., Hepburn, A., Kirk, J. A., & Kekenes-Huskey, P. M. (2025). GSK3β-driven phosphorylation of ABLIM1 regulates its interactions with titin in cardiac muscle. *Journal of General Physiology, 157*(5), e202413737. [https://doi.org/10.1085/jgp.202413737](https://doi.org/10.1085/jgp.202413737)

## 🙏 Acknowledgements

I'd like to thank Dr. Kekenes-Huskey and Alec Loftus at the Stritch School of Medicine for their guidance and support throughout this project.
