# Sampling Algorithms for High Average-Efficiency Itemsets under Length Constraints

**Cite:** *"Ashraf, M., Amin, E., Diop, L., Plantevit, M., Abdelkader, T., Ismail, R. (2026). Sampling algorithms for high average-efficiency itemsets under length constraints. Knowledge-Based Systems."* (to appear — the published reference and the PDF will be added here once available)

High-efficiency pattern mining algorithms discover itemsets whose return is high relative to the investment they require, which makes them more informative than high utility itemsets whenever the action triggered by a pattern has a measurable cost. Existing algorithms enumerate the complete set of patterns above a minimum efficiency threshold, so their output is hypersensitive to that threshold and tends to focus on a few overlapping groups. We propose two algorithms, HAES-UB and HAES-MH, that sample itemsets so that each itemset is drawn with a probability proportional to its average efficiency in the database, under length constraints that discard long, rare itemsets composed of low-investment items. HAES-UB is a rejection sampler built on the Efficiency Bound Matrix and three upper bounds on average efficiency; HAES-MH is an independence Metropolis-Hastings chain that reuses HAISampler as its proposal and corrects it with an investment ratio test. Neither algorithm requires a minimum efficiency threshold. Experiments show that both extract thousands of high-average-efficiency patterns in a few seconds from different databases.

## Content

| File | Description |
|---|---|
| `HAES-UB.py` | Rejection sampler. The bound-in-use is chosen with `v` = 1, 2, or 3 (`eub1`, `eub2`, `eub3`). |
| `HAES-MH.py` | Independence Metropolis-Hastings sampler. The HAISampler proposal is embedded in the file. |
| `DatasetsHEP/` | Databases and their investment values. |
| `XP_Output/` | Sampled patterns, written by the two algorithms. |

## Input format

Each database is a text file in SPMF format, one transaction per line, with three fields separated by `:`

```
items: transaction utility: utility of each item
1 3 7: 26: 5 8 13
```

Each database is accompanied by a file `<name>_investmentValues.txt` giving the total investment of every item, one item per line:

```
1 3731
3 455
7 10283
```

The two files of a database `foodmart` are therefore `DatasetsHEP/foodmart.txt` and `DatasetsHEP/foodmart_investmentValues.txt`.

## Usage

Set the parameters at the top of the `main` function and run the file:

```python
datasets = ["foodmart"]   # databases to sample, taken from DatasetsHEP
N = 10000                 # the desired sample size
m = 1                     # minimum length constraint
M = 5                     # maximum length constraint
v = 3                     # HAES-UB only: the upper bound in use, 1, 2 or 3
```

```
python HAES-UB.py
python HAES-MH.py
```

Both print the preprocessing time, the sampling time, the acceptance rate and the number of distinct patterns, and write the sample to `XP_Output/`. Each output file lists one pattern per line with its frequency in the sample:

```
Pattern	&	Frequency in the sample
534	36
1096 1342 321 327 664	1
```

## Requirements

Python 3, standard library only.

## License

This project is released under the GNU General Public License v3.0. `HAES-MH.py` embeds the weighting and drawing routines of HAISampler (Diop, PAKDD 2022, GPL-3.0), reused here as the proposal distribution:

> Diop, L. (2022). High Average-Utility Itemset Sampling Under Length Constraints. In: Advances in Knowledge Discovery and Data Mining. PAKDD 2022. Lecture Notes in Computer Science, vol 13281. Springer, Cham. https://doi.org/10.1007/978-3-031-05936-0_11
