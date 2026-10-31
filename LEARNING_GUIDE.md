# Learning guide - dead-simple explanations

This guide explains every idea in this project as simply as possible. Read it top to
bottom once, then open the code file named in each section and match the words to the lines.

---

## 1. What is "alignment"?

Two DNA pieces that came from the same ancestor look *almost* the same:

```text
GATTACA
GCATGCU
```

Alignment means lining them up, adding gaps (`-`) where letters were inserted or lost,
so similar letters sit on top of each other:

```text
G-ATTACA
| | |.|.        | = same letter, . = different letter, space = gap
GCA-TGCU
```

We give points to each column (**scoring**): for example +1 for same, -1 for different,
-1 for a gap. The alignment with the **highest total** is the best one.

The problem: there are *astronomically* many ways to place gaps. Trying them all is
impossible. **Dynamic programming** finds the best one quickly.

---

## 2. Dynamic programming (DP) in one minute

DP = "solve big problems by reusing answers to small problems, written in a table."

Think of climbing stairs: the number of ways to reach step 10 = ways to reach step 9 +
ways to reach step 8. If you write each answer in a list, you never compute it twice.

For alignment, the small problem is: *"what is the best score for aligning the first `i`
letters of A with the first `j` letters of B?"* Call it `S[i][j]`.

The last column of that alignment can only be one of three things:

1. `A[i]` on top of `B[j]` -> `S[i-1][j-1] + score(A[i], B[j])`  (diagonal arrow)
2. `A[i]` on top of a gap  -> `S[i-1][j] + gap`                     (up arrow)
3. a gap on top of `B[j]`  -> `S[i][j-1] + gap`                     (left arrow)

So `S[i][j] = max(those three)`. We fill the table row by row and remember which arrow
won. That is the whole algorithm. Code: `_linear()` in `src/dnaalign/align.py`.

---

## 3. The matrices (tables), drawn out

A = `GAT`, B = `GT`, match +1, mismatch -1, gap -1:

```text
          ""    G     T
    ""     0   -1    -2        first row/column = "all gaps so far"
    G     -1    1     0
    A     -2    0     0
    T     -3   -1     1   <- bottom-right = best global score = 1
```

**Traceback**: start at the bottom-right cell and follow the remembered arrows back to the
top-left. Each arrow gives one column of the alignment, written backwards:

```text
GAT
| |
G-T      score = 1 + (-1) + 1 = 1
```

### Global vs local

- **Needleman-Wunsch (global)**: align *everything*. Answer = bottom-right cell.
- **Smith-Waterman (local)**: find the best *matching piece*. Two changes:
  1. a cell may never go below 0 (0 means "start fresh here"),
  2. the answer is the **largest cell anywhere**, and traceback stops at a 0.

Use local when one sequence is a small piece of a big one (a read inside a genome).

### Affine gaps (Gotoh) - why 3 tables?

In real biology one long gap of 3 letters is much more likely than 3 separate gaps.
So we charge a big **open** cost and a small **extend** cost:
gap of length L = `open + (L-1) * extend` (e.g. -5, -1, -1 = -7 for length 3).

To know whether a gap is new or continuing, we keep 3 tables:

- `M` - this column is letter vs letter
- `X` - this column is a letter of A vs a gap (we were "inside a gap in B")
- `Y` - this column is a gap vs a letter of B

Moving inside `X` costs `extend`; jumping into `X` from `M` or `Y` costs `open`.
Code: `_affine()` in `align.py`.

### Scoring matrices

- DNA: simple match/mismatch (`Scoring(match=2, mismatch=-4, ...)`).
- Proteins: some amino-acid swaps are common and harmless (e.g. I <-> V), others are rare.
  **BLOSUM62** is a 24x24 table of log-odds scores learned from real protein families.
  Positive = "seen together often", negative = "rarely". File: `src/dnaalign/data/BLOSUM62.txt`.

---

## 4. Big-O: why alignment is O(n·m)

Big-O answers: *"if the input gets bigger, how much slower does it get?"*

- The table has `(n+1) x (m+1)` cells.
- Each cell does a fixed amount of work (3 additions, 1 max).
- So time = constant x n x m -> **O(n·m)**. Memory is also O(n·m) because we keep arrows.

Our measured proof (`figures/dp_scaling.png`): length 200 took 0.031 s, length 400 took
0.137 s, length 800 took 0.587 s. Doubling the length -> about 4x the time (2 x 2).

Why that matters: a 150-letter read against a 30,000-letter genome = 4.5 million cells per
read. Against the human genome (3 billion letters) = 450 billion cells per read. Too slow!
That is why real mappers use an **index** first.

| Big-O | Name | Example |
|---|---|---|
| O(1) | constant | dictionary lookup |
| O(log n) | logarithmic | binary search |
| O(n) | linear | building the k-mer index |
| O(n log n) | | sorting |
| O(n·m) / O(n^2) | quadratic | DP alignment |

---

## 5. Hashing and k-mers

A **k-mer** is a word of `k` letters. The 3-mers of `ACGTACG` are `ACG CGT GTA TAC ACG`.

A **hash map** (Python `dict`) stores `key -> value` and finds a key in O(1) on average:
it turns the key into a number (the *hash*) and jumps straight to that slot, instead of
searching one by one.

Our **k-mer index** is `dict[k-mer] -> list of positions`:

```text
ACG -> [0, 4]
CGT -> [1]
```

Building it: one pass over the genome -> O(L). File: `src/dnaalign/kmer_index.py`.

### Seed-and-extend (how the mapper works) - `mapper.py`

1. **Seed**: cut the read into k-mers, look each one up (fast).
2. **Vote**: a hit "read position 10 = genome position 5010" means the read starts at
   genome position `5010 - 10 = 5000` (its *diagonal*). Correct hits all vote for the
   same diagonal; random hits scatter.
3. **Extend**: run Smith-Waterman only in a thin strip (**band**) around the winning
   diagonal: ~150 x 33 cells instead of 150 x 30,000 -> about 900x less work.
4. Try the reverse complement too, because DNA has two strands.

Choosing `k`: big `k` = fewer random hits (fast) but one sequencing error breaks a seed
(less sensitive). Our benchmark shows it: at 15% errors, k=15 placed 90.5% of reads,
k=11 placed 99.5% but was ~3.3x slower.

---

## 6. Likely interview questions (with short answers)

**Q1. Explain Needleman-Wunsch.**
DP table where `S[i][j]` = best score aligning first `i` letters of A with first `j` of B.
Each cell = max of diagonal + match/mismatch, up + gap, left + gap. Answer at the
bottom-right; traceback gives the alignment. O(n·m) time and memory.

**Q2. Difference between global and local alignment?**
Global aligns full sequences (answer at the corner). Local (Smith-Waterman) clamps
scores at 0 and takes the maximum anywhere, so it finds the best matching region.

**Q3. Why affine gaps and why three matrices?**
Biology: one long gap is more likely than many short ones. To charge open vs extend
differently you must remember whether you are already in a gap, so you keep one table per
state (match, gap in A, gap in B). Still O(n·m).

**Q4. Can you reduce memory?**
If you only need the score, keep just two rows -> O(m) memory. To also get the alignment
in linear memory use Hirschberg's divide-and-conquer.

**Q5. How do you speed up alignment against a whole genome?**
Index the genome (k-mer hash, suffix array, or FM-index), find exact seeds, then align
only near the seeds, in a band. That is seed-and-extend (BLAST, BWA, minimap2).

**Q6. What is a hash map and what is its lookup cost?**
An array of buckets addressed by `hash(key)`. Average O(1) insert/lookup; worst case O(n)
if many keys collide.

**Q7. How do you pick k?**
Large enough that random matches are rare (4^k much bigger than the genome length:
4^15 ~ 1 billion vs 30 kb), small enough that an error-free k-mer still exists in a noisy
read. Trade-off: speed vs sensitivity (see results table).

**Q8. How did you know your code was correct?**
I compared scores with Biopython's PairwiseAligner on 3,000 random cases (3,000/3,000
identical), re-scored every traceback, wrote 60+ pytest tests, and ran CI on Python 3.10-3.13.

**Q9. How did you measure the mapper?**
Simulated reads from the real SARS-CoV-2 genome with known start positions and 0-15%
errors, then counted reads placed on the right strand within 10 bp.

**Q10. What would you do next?**
Rewrite the inner loop in C/Numba or with SIMD, add chaining of multiple seeds, output
SAM, support multiple chromosomes, and test on real Illumina reads.

**Coding warm-ups to practise** (LeetCode): Edit Distance (#72), Longest Common
Subsequence (#1143), Two Sum (#1, hash map), Repeated DNA Sequences (#187, k-mers + hashing).
