# TreesAndModels — Genomics Assignment

This folder contains one assignment with two parts:

- **Part A — Sequence modeling & classification:** first-order Markov modeling, coding vs. neutral scoring using 64×64 probability matrices, and ROC analysis.
- **Part B — UPGMA tree building:** hierarchical clustering from scratch.

The code is self-contained so it can live neatly inside a larger `genomics` repository.

---

## Contents

- README.md: This file;
- Ancestor.fa: Ancestor/reference FASTA (if required by your scripts);
- Spacii.fa: FASTA set 1;
- Spacii_2100.fa: FASTA set 2;
- codingModel.tab: 64×64 probability matrix (coding model);
- noncodingModel.tab: 64×64 probability matrix (neutral/noncoding model);
- classifierspacii2100.py: first-order Markov classifier (Part A);
- codingModelClassifier.py: coding vs. neutral scoring with 64×64 matrices (Part A);
- ROC.py: ROC/AUC utilities/plotting (Part A);
- rates.py: Shared metrics helpers (TP/FP/FN/TN, accuracy, etc.);
- UPGMA.py: UPGMA implementation (Part B).

> **Note:** Some scripts may currently expect files to be alongside them (no CLI flags). If a path error occurs, open the script and adjust the file names/paths at the top, or add simple `argparse` flags later.

---

## Environment

- **Python:** 3.9+ (3.10+ recommended)
- **Packages:** `numpy`, `pandas`, `matplotlib`, `biopython`

Set up a quick virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install numpy pandas matplotlib biopython
```

## Part A — Sequence Modeling & Classification

### What Part A does

- **First-order Markov classifier (`classifierspacii2100.py`)**
  - Learns transition probabilities `P(x_i | x_{i-1})` from training sequences.
  - Scores test sequences by total log-likelihood under each class model and predicts the class with the higher likelihood.

- **Coding vs. neutral scoring (`codingModelClassifier.py`)**
  - Uses 64×64 codon matrices (`codingModel.tab`, `noncodingModel.tab`) to score sequences under coding vs. neutral assumptions.
  - Produces per-sequence log-likelihoods or likelihood ratios that you can threshold to classify.

- **ROC/AUC & metrics (`ROC.py`, `rates.py`)**
  - Given true labels and continuous scores, computes TPR/FPR across thresholds and AUC; optionally plots an ROC curve.

### Run examples

```bash
# 1) Train/apply the first-order Markov classifier
python classifierspacii2100.py

# 2) Score under coding vs neutral models (64×64)
python codingModelClassifier.py

# 3) Compute/plot ROC (adapt to how your scripts pass/expect scores)
python ROC.py
```

### Outputs to expect

- **Scores/Predictions:** CSV or text table with per-sequence scores (e.g., log-likelihoods or ratios).
- **ROC/AUC:** CSV with TPR/FPR per threshold and/or an ROC plot (PNG).
- **Threshold choice:** Choose a cutoff (e.g., maximizing AUC or targeting ~95% specificity) and note FP/FN trade-offs.

## Part B — UPGMA Tree Building

### What Part B does

- **Inputs are hard-coded** at the top:
  - `distanceMatrix`: a symmetric matrix with zeros on the diagonal.
  - `speciesList`: a list of labels (same order as the matrix).
- **Finds the closest pair** using `findSmallest(dM)`, scanning the upper triangle for the **smallest positive** off-diagonal entry.
  - **Tie-breaking is random**: if multiple pairs share the minimum distance, one is chosen at random.
- **Merges clusters** via `updateMatrix(dM, row, col)` using a **simple (unweighted) average** of distances from the two merged clusters to every other cluster.
  - *Note:* This is a plain arithmetic mean, **not size-weighted** by cluster sizes.
- **Updates labels** via `updateSpecies(sp, r, c, branch_length)` by replacing the two merged labels with a string of the form:
  - `(<left>,<right>):<branch_length>`
  - Here, `branch_length` is set to the **pairwise distance** between the merged clusters at the time of merging.
- **Prints intermediate state** at each iteration:
  - `"Updated Distance Matrix:"` followed by the new matrix.
  - `"Updated Species List:"` showing the progressively nested label strings.


### Run example

```bash
python UPGMA.py
```

### Outputs to expect

- A series of iterations, each printing:
  - **Updated Distance Matrix:** the reduced/updated matrix after merging the closest pair.
  - **Updated Species List:** the list of labels with the merged pair replaced by a parenthetical string (e.g., `"(T_Pain,G_Unit):2"`).
