# Project Report - Alignment

## Baseline

### Design Experience

I worked with Brandon Monson and Porter Schollenberger to design an unrestricted global alignment using Needleman–Wunsch. We planned a full (n+1)×(m+1) dynamic programming table with separate backpointers, simple linear gap penalties (no affine open), and a traceback from (n,m) to build the aligned strings. We chose clarity over micro-optimizations to make correctness and debugging straightforward.

### Theoretical Analysis - Unrestricted Alignment

#### Time

```py
n, m = len(seq1), len(seq2)

# allocate DP and backpointer tables
editCostTable = [[0]*(m+1) for _ in range(n+1)]                                   #            O(nm): allocate and zero-initialize
backPointerTable = [[""]*(m+1) for _ in range(n+1)]                                #            O(nm): allocate backpointers

# initialize first column / first row (global alignment borders)
for i in range(1, n + 1):                                                          #            O(n)
    editCostTable[i][0] = editCostTable[i-1][0] + indel_penalty
    backPointerTable[i][0] = "U"
for j in range(1, m + 1):                                                          #            O(m)
    editCostTable[0][j] = editCostTable[0][j - 1] + indel_penalty
    backPointerTable[0][j] = "L"

# fill DP matrix
for i in range(1, n + 1):                                                          #            O(n)
    for j in range(1, m + 1):                                                      #            O(m) → nested ⇒ O(nm)
        diagonal = editCostTable[i-1][j-1] + (match_award if seq1[i-1]==seq2[j-1]
                                              else sub_penalty)                    #            O(1)
        left     = editCostTable[i][j-1] + indel_penalty                           #            O(1)
        up       = editCostTable[i-1][j] + indel_penalty                           #            O(1)
        # choose min and record backpointer
        cost, _, direction = min([(diagonal,0,"D"), (left,1,"L"), (up,2,"U")])     #            O(1)
        editCostTable[i][j] = cost
        backPointerTable[i][j] = direction

# traceback to build alignments
i, j = n, m
aligned1, aligned2 = [], []
while i>0 or j>0:                                                                  #            O(n+m): one step per move
    d = backPointerTable[i][j]
    if d == 'D': aligned1.append(seq1[i-1]); aligned2.append(seq2[j-1]); i-=1; j-=1
    elif d == 'L': aligned1.append(gap);        aligned2.append(seq2[j-1]); j-=1
    else:        aligned1.append(seq1[i-1]);    aligned2.append(gap);        i-=1
```

The DP fill dominates with **O(nm)** time; initialization adds **O(n+m)** and traceback adds **O(n+m)**, which are lower-order terms. Overall runtime is **O(nm)**.

#### Space

```py
# main tables
editCostTable     = [[0]*(m+1) for _ in range(n+1)]                                #            O(nm): numeric DP table
backPointerTable  = [[""]*(m+1) for _ in range(n+1)]                                #            O(nm): backpointers

# output strings built by traceback
aligned1_list: list[str] = []                                                      #            O(n+m): stores final alignment 1
aligned2_list: list[str] = []                                                      #            O(n+m): stores final alignment 2
```

Ignoring the input sequences, auxiliary memory is dominated by the DP and backpointer tables for **O(nm)**, plus **O(n+m)** for the outputs; total auxiliary space is **O(nm)**. Including the input (and outputs), overall storage remains **O(nm)**.

### Empirical Data - Unrestricted Alignment

| N    | Time (sec) |
| ---- | ---------- |
| 500  | 0.076      |
| 1000 | 0.318      |
| 1500 | 0.735      |
| 2000 | 1.304      |
| 2500 | 2.058      |
| 3000 | 2.96       |

### Comparison of Theoretical and Empirical Results - Unrestricted Alignment

- Theoretical order of growth: **O(nm)**
- Empirical order of growth (if different from theoretical):

![](fill-me-in.png)

_Fill me in_

## Core

### Design Experience

_Fill me in_

### Theoretical Analysis - Banded Alignment

#### Time

_Fill me in_

#### Space

_Fill me in_

### Empirical Data - Banded Alignment

| N     | time (ms) |
| ----- | --------- |
| 100   |           |
| 1000  |           |
| 5000  |           |
| 10000 |           |
| 15000 |           |
| 20000 |           |
| 25000 |           |
| 30000 |           |

### Comparison of Theoretical and Empirical Results - Banded Alignment

- Theoretical order of growth:
- Empirical order of growth (if different from theoretical):

![](fill-me-in.png)

_Fill me in_

### Relative Performance Of Unrestricted Alignment versus Banded Alignment

_Fill me in_

## Stretch 1

### Design Experience

_Fill me in_

### Code

```python
# Fill me in
```

### Alignment Scores

_Fill me in_

## Stretch 2

### Design Experience

_Fill me in_

### Empirical Data - Using Affine Penalties

| N    | time (ms) |
| ---- | --------- |
| 500  |           |
| 1000 |           |
| 1500 |           |
| 2000 |           |
| 2500 |           |
| 3000 |           |

### Empirical Outcome Comparisons

_Fill me in_

### Alignment Outcome Comparisons

##### Sequences and Alignments

_Fill me in_

##### Chosen Parameters and Better Alignments Discussion

_Fill me in_

## Project Review

_Fill me in_
