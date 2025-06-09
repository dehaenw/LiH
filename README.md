# LiH
Ligand hallucinations in cofolding methods. Methods to quantify common ligand hallucinations in cofolding methods. In particular, methods to detect non-planar rings that should be planar and vice versa. 

# Introduction
On common protein-ligand pose prediction benchmarks, cofolding methods such as alphafold 3 and boltz outperform or match molecular docking, the most commonly used traditional pose prediction method.

Nonetheless, some unusual issues plague these new methods. The authors of AlphaFold3 already point out that asymmetric carbons are frequently inverted (4.4%) in the output of their method, but otherwise they manage to pass standard validity checks (PoseBusters). Surprisingly, we have found several examples where the input molecules (given as SMILES) did not match the output molecules (given as CIF). Some common problem motifs are: cyclohexanes (aromatized to benzenes), tetrahydrofurans (aromatized to furan), allenes, alkynes and nitriles (all 3 with non-linear structures). These can be considered an example of the “hallucination” phenomenon which is well known in other generative AI methods such as LLMs. These ligand hallucinations are also interesting because they subvert pose validation methods such as RMSD based checks compared to a reference.

We propose a method to systematically identify instances of cofolding ligand hallucinations by generating ligand structures, inferring their connectivity (using the input structure as a template) and comparing it to the input connectivity. We use this to compile a list of problematic functionalities which can be used to compare the extent of ligand hallucinations across cofolding methods.

# Example use
The jupyter notebook `Calculate_Example.ipynb` contains an example calculation of the hallucination rate for ETKDG+MMFF94(the reference) versus Vina (precalculated and provided as sdf) versus just ETKDG versus Boltz (precalculated and provided as sdf). This gives the following output, highlighting the higher hallucination rate of Boltz compared to traditional methods.

# Working repo
This is a preview version. Expect a more extended repo later, this repo just provides tha current version of the cofolding tricky ligand set as well as methods for conformation assessment. 

# LiH set
The set `lih_set.py` counts 100 molecules and is divided in 4 groups:
- Collinear set which contains atom sets that should be collinear (e.g. allenes, alkynes, ...)
- Non-planar molecules such as cyclobutanes,cyclopentanes and cyclohexanes and their heterocyclic analogs
- Planar molecules such as aromatic rings, in particular those decorated with hydroxy groups (which can sometimes lead to wrong sugar/inositol like predict conformations)
- Drugs and druglike molecules from drugbank
  
changes of the content and size of this set are still expected.
