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
- Empirical order of growth (if different from theoretical): N/A

![img](./_analysis/align_empirical_unbounded.svg)

The empirical and observed match up almost perfectly.

## Core

### Design Experience

I worked with Brandon Monson and Porter Schollenberger to design a banded Needleman–Wunsch. We restrict DP to cells with `|i−j| ≤ d`, compute row-specific column bounds via `j_bounds`, and keep rolling arrays for costs and pointers while saving one backpointer slice per row for traceback. By hand, we walked a small example to verify band entry/exit and diagonal/left/up transitions. We chose arrays (lists) over dicts for cache-friendly scans across the band; the tradeoff is storing `O(n·d)` backpointers to enable simple traceback.

### Theoretical Analysis - Banded Alignment

#### Time

```py
n, m = len(seq1), len(seq2)
d = max(n, m) if banded_width < 0 else banded_width                               #            sets band half-width d

if banded_width >= 0 and abs(n - m) > banded_width:                               #            O(1): early reject when paths can't fit band
    return inf, None, None

# rolling rows over the band (width ≤ 2d+1)
prev_cost = [0] * (2 * d + 1)                                                     #            O(d)
curr_cost = [0] * (2 * d + 1)                                                     #            O(d)
prev_ptr  = [""] * (2 * d + 1)                                                    #            O(d)
curr_ptr  = [""] * (2 * d + 1)                                                    #            O(d)

# initialize i = 0 row within band columns
lo, hi = j_bounds(0, m, d)                                                        #            O(1)
for j_offset, j in enumerate(range(lo, hi + 1)):                                   #            O(d): init banded first row
    prev_cost[j_offset] = j * indel_penalty
    prev_ptr[j_offset]  = "L" if j > 0 else ""

backpointer_rows = [prev_ptr.copy()]                                              #            O(d)

# fill DP within band
for i in range(1, n + 1):                                                         #            n iterations
    lo, hi = j_bounds(i, m, d)                                                    #            O(1)
    prev_lo, prev_hi = j_bounds(i - 1, m, d)                                      #            O(1)
    for j_offset, j in enumerate(range(lo, hi + 1)):                              #            O(d) cells per row → O(n·d)
        candidates = []
        if prev_lo <= j-1 <= prev_hi:                                             #            O(1)
            diag_idx = (j - 1) - prev_lo
            diagonal = prev_cost[diag_idx] + (match_award if seq1[i-1]==seq2[j-1]
                                              else sub_penalty)                   #            O(1)
            candidates.append((diagonal, 0, 'D'))
        if j - 1 >= lo:                                                            #            O(1)
            left = curr_cost[j_offset - 1] + indel_penalty
            candidates.append((left, 1, "L"))
        if prev_lo <= j <= prev_hi:                                                #            O(1)
            up_idx = j - prev_lo
            up = prev_cost[up_idx] + indel_penalty
            candidates.append((up, 2, "U"))

        min_cost, _, direction = min(candidates)                                   #            O(1)
        curr_cost[j_offset] = min_cost
        curr_ptr[j_offset]  = direction

    backpointer_rows.append(curr_ptr.copy())                                       #            O(d) per row
    prev_cost, curr_cost = curr_cost, [0] * (2 * d + 1)                            #            O(d)
    prev_ptr,  curr_ptr  = curr_ptr,  [""] * (2 * d + 1)                           #            O(d)

# traceback walks at most n+m steps, staying inside band
i, j = n, m
aligned1_list, aligned2_list = [], []
while (i > 0 or j > 0):                                                            #            O(n+m)
    lo, hi = j_bounds(i, m, d)                                                     #            O(1)
    j_offset = j - lo                                                              #            O(1)
    dch = backpointer_rows[i][j_offset]                                            #            O(1)
    if dch == 'D': i -= 1; j -= 1
    elif dch == 'L': j -= 1
    else: i -= 1
```

Each of the `n` rows processes at most `~2d+1` cells, so the DP fill is **O(n·d)** (symmetrically **O(m·d)**; for `n≈m`, **O(n·d)**). Initialization is **O(d)** and traceback **O(n+m)**, which are lower-order compared to **O(n·d)**.

#### Space

```py
# rolling band rows
prev_cost = [0] * (2 * d + 1)                                                     #            O(d)
curr_cost = [0] * (2 * d + 1)                                                     #            O(d)
prev_ptr  = [""] * (2 * d + 1)                                                    #            O(d)
curr_ptr  = [""] * (2 * d + 1)                                                    #            O(d)

# store backpointers for each i to enable traceback
backpointer_rows: list[list[str]] = [prev_ptr.copy()]                             #            O(d) initially
# ...
backpointer_rows.append(curr_ptr.copy())  # per row                               #            O(n·d) total across all rows

# output buffers
aligned1_list: list[str] = []                                                     #            O(n+m)
aligned2_list: list[str] = []                                                     #            O(n+m)
```

Ignoring the input sequences, auxiliary memory is dominated by stored backpointer slices at **O(n·d)**, plus rolling rows **O(d)** and outputs **O(n+m)**, so total auxiliary space is **O(n·d)**. Including the input and outputs, the overall storage remains **O(n·d)**.

### Empirical Data - Banded Alignment with band_width = 5

| N    | Time (sec) |
| ---- | ---------- |
| 500  | 0.002      |
| 1000 | 0.005      |
| 1500 | 0.008      |
| 2000 | 0.01       |
| 2500 | 0.013      |
| 3000 | 0.02       |

### Comparison of Theoretical and Empirical Results - Banded Alignment

- Theoretical order of growth: **O(nd)** where d is the band width
- Empirical order of growth (if different from theoretical):

![img](./_analysis/align_empirical_bounded.svg)

The empirical lines up pretty well with the expected, however there are some random outliers at the end, and I'm not sure what those are or why they randomly spike so much higher for some reason.

### Relative Performance Of Unrestricted Alignment versus Banded Alignment

Just based on the empirical data tables with a band of 5, it is apparent how much faster the banded width algorithm really is. The graphs side by side show that one experiences exponential growth while the other is linear because the band width is constant.

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
