---
title: Strassen–Winograd multiplication
summary: The 7-multiplication scheme FLINT uses for large matrix products over $\Z/n\Z$, with Bodrato's operation sequence, a schedule with two temporaries, and odd dimensions handled by virtual zero padding.
complexity: $O(d^{\log_2 7}) \approx O(d^{2.807})$ operations for $d \times d$ matrices; one level costs 7 products of half size and 15 additions of quarter-size blocks.
order: 15
---

Strassen showed that two $2 \times 2$ matrices over any ring can be
multiplied with 7 multiplications instead of 8 [@Str1969]. Since the
entries can themselves be matrices, applying the scheme recursively
multiplies $d \times d$ matrices in $O(d^{\log_2 7})$ ring operations.
Winograd's variant [@Win1971] needs 15 block additions instead of 18.
The version described below is that of `nmod_mat_mul_strassen`, which
is also used, in byte or bit-packed form, by `nmod_mat_mul_u8`.

## Generic version {#generic}

`gr_mat_mul_strassen` performs one level of the Winograd form of
Strassen's algorithm over any ring, in Bodrato's order of operations,
computing the seven block products with `gr_mat_mul` and handling odd
dimensions by peeling off the last row or column. Recursion therefore
happens only through the ring's matrix multiplication method: a ring
installs a method that calls `gr_mat_mul_strassen` above a cutoff and a
base-case product below it, as the context `nmod32` does on top of
the 32-bit SIMD kernel of `nmod_mat`. The scheme needs no commutativity, so it also
applies to matrices over noncommutative rings and to block matrices. The
same seven-product formulas multiply the $2 \times 2$ matrices of
polynomials in the half-gcd algorithm of `gr_poly`
([](algorithms/polynomials/half-gcd)) once their entries have length
20 or more.

## The formulas {#formulas}

Split $A$ ($m \times k$), $B$ ($k \times n$) and $C = AB$ into $2 \times 2$
blocks $A_{ij}$, $B_{ij}$, $C_{ij}$ of compatible sizes. FLINT uses the
sequence of Bodrato [@Bodrato2010], which is equivalent to Winograd's
and was designed so that the same operands serve for squaring. With

$$
\begin{aligned}
S &= A_{12} + A_{22} - A_{21}, & T &= B_{12} + B_{22} - B_{21}, \\
P_1 &= (A_{12} + A_{22})(B_{12} + B_{22}), & P_2 &= (A_{22} - A_{21})(B_{22} - B_{21}), \\
P_3 &= S \, T, & P_4 &= (S - A_{11}) \, B_{12}, \\
P_5 &= A_{12} B_{21}, & P_6 &= A_{21} (T - B_{11}), \\
P_7 &= A_{11} B_{11},
\end{aligned}
$$

the blocks of the product are

$$
\begin{aligned}
C_{11} &= P_5 + P_7, & C_{12} &= P_3 + P_5 - P_4 - P_2, \\
C_{21} &= P_1 - P_3 - P_5 - P_6, & C_{22} &= P_1 + P_2 - P_3 - P_5.
\end{aligned}
$$

::: theorem {#thm-correct}
These formulas give $C = AB$ over any ring, commutative or not.
:::

::: proof
$C_{11} = A_{11} B_{11} + A_{12} B_{21}$ is immediate. Let
$U = A_{12} + A_{22}$ and $V = B_{12} + B_{22}$, so $S = U - A_{21}$ and
$T = V - B_{21}$. Then

$$
P_1 - P_3 = UV - (U - A_{21})(V - B_{21}) = U B_{21} + A_{21} V - A_{21} B_{21}.
$$

Subtracting $P_6 = A_{21} V - A_{21} B_{21} - A_{21} B_{11}$ and
$P_5 = A_{12} B_{21}$ leaves $A_{22} B_{21} + A_{21} B_{11} = C_{21}$.
Adding $P_2 - P_5 = (A_{22} - A_{21})(B_{22} - B_{21}) - A_{12} B_{21}$
to $P_1 - P_3$ instead gives, after cancellation,
$A_{21} B_{12} + A_{22} B_{22} = C_{22}$. Finally
$P_3 - P_4 = S T - S B_{12} + A_{11} B_{12} = S (B_{22} - B_{21}) + A_{11} B_{12}$
and $S (B_{22} - B_{21}) = A_{12} (B_{22} - B_{21}) + P_2$, so
$P_3 - P_4 - P_2 + P_5 = A_{11} B_{12} + A_{12} B_{22} = C_{12}$. No
step exchanges the two factors of a product, so commutativity is not
used.
:::

The operands of the products need 8 additions ($U$, $V$, $A_{22} - A_{21}$,
$B_{22} - B_{21}$, $S$, $T$, $S - A_{11}$, $T - B_{11}$) and the results
7 more, 15 in all.

## Schedule and memory {#schedule}

The implementation follows the memory-efficient schedule of Dumas,
Pernet and Zhou [@DPZ2009]: the products are written into the blocks of
$C$ as they are computed, and only two temporaries are needed, $X_1$ of
the size of $A_{11}$ (or $C_{11}$, whichever is wider) and $X_2$ of the
size of $B_{11}$. With $X_1, X_2$ and the four blocks of $C$ as storage:

::: algorithm Strassen–Winograd, one level {#alg-sw}
1. $X_1 \gets A_{22} + A_{12}$, $X_2 \gets B_{22} + B_{12}$,
   $C_{21} \gets X_1 X_2$ ($P_1$).
2. $X_1 \gets A_{22} - A_{21}$, $X_2 \gets B_{22} - B_{21}$,
   $C_{22} \gets X_1 X_2$ ($P_2$).
3. $X_1 \gets X_1 + A_{12}$, $X_2 \gets X_2 + B_{12}$,
   $C_{11} \gets X_1 X_2$ ($P_3$).
4. $X_1 \gets X_1 - A_{11}$, $C_{12} \gets X_1 B_{12}$ ($P_4$).
5. $X_1 \gets A_{12} B_{21}$ ($P_5$), $C_{11} \gets C_{11} + X_1$.
6. $C_{12} \gets C_{12} + C_{22}$, $C_{12} \gets C_{11} - C_{12}$ (final).
7. $C_{11} \gets C_{21} - C_{11}$ ($P_1 - P_3 - P_5$).
8. $X_2 \gets X_2 - B_{11}$, $C_{21} \gets A_{21} X_2$ ($P_6$).
9. $C_{21} \gets C_{11} - C_{21}$, $C_{22} \gets C_{22} + C_{11}$ (final).
10. $C_{11} \gets A_{11} B_{11}$ ($P_7$), $C_{11} \gets X_1 + C_{11}$ (final).
:::

The 7 products are computed by `nmod_mat_mul`, so they use the fast
kernels or recurse into Strassen–Winograd again, and $C$ must not be
aliased with $A$ or $B$. If any dimension is at most 4,
`nmod_mat_mul_strassen` calls `nmod_mat_mul` directly.

## Odd dimensions {#odd}

FLINT splits a dimension $d$ into $h = \lceil d/2 \rceil$ for the top or
left blocks and $d - h \in \{h - 1, h\}$ for the bottom or right blocks.
Conceptually, the short blocks are padded with one zero row or column, so
that all blocks have the same size and [](#thm-correct) applies.
Padding $A$ with a zero row and a zero column (and $B$ accordingly)
pads $AB$ with zeros, so the padded product contains $C$ as its top left
block. Nothing is actually padded:

- the temporaries $X_1$, $X_2$ have the full padded size, and the sums of
  blocks that form them copy the extra row or column of the larger
  operand and write zeros beyond;
- a product with a short block of $A$ or $B$ as an operand uses the real
  dimensions of that block;
- the blocks $C_{12}$, $C_{21}$ and $C_{22}$ of the real $C$ lack the
  padding, while the schedule stores intermediate products in them.

::: proposition {#prop-odd}
With the schedule of [](#alg-sw), the entries missing from $C_{12}$,
$C_{21}$, $C_{22}$ are never needed for an entry of $C$.
:::

::: proof
If $n$ is odd, $C_{12}$ and $C_{22}$ lack their last column. The final
$C_{11} = P_5 + P_7$ and $C_{21} = P_1 - P_3 - P_5 - P_6$ do not read
$C_{12}$ or $C_{22}$, and the final $C_{12}$ and $C_{22}$ only need
their own columns, so the missing column is not used. If $m$ is odd,
$C_{21}$ and $C_{22}$ lack their last row. The last rows of the values
held there ($P_1$, $P_2$, $P_6$ and the partial sums of step 9) only feed
the last rows of the final $C_{21}$ and $C_{22}$, which do not exist,
except for $P_2$, which is added to $C_{12}$ in step 6. But
$P_2 = (A_{22} - A_{21})(B_{22} - B_{21})$ has a zero last row, because
the padded $A_{21}$ and $A_{22}$ both have a zero last row. If $k$ is
odd, only the temporaries and the inner dimension of $P_5 = A_{12} B_{21}$
are affected, and both are handled exactly as stated above.
:::

So odd dimensions cost nothing beyond the padded additions, which touch
the same entries as ordinary ones. Earlier versions of FLINT peeled off
the last row or column and corrected the result with matrix-vector
products, a technique known as dynamic peeling [@HLJJTG1996].

## Complexity {#complexity}

Let $T(d)$ be the cost for $d \times d$ matrices when the recursion goes
down to dimension $d_0$. Then
$T(d) = 7 T(\lceil d/2 \rceil) + 15 \lceil d/2 \rceil^2$, so
$T(d) = O(d^{\log_2 7})$. One level replaces 8 products by 7, saving one
eighth of the multiplication work, at the cost of 15 additions of
$(d/2)^2$ entries, each of which reads and writes memory. The level pays
off when a product of dimension $d/2$ costs more than about 15 additions
of that size. With the classical algorithm, which forms one dot product
per entry, this happens around dimension 200 to 400 in FLINT. With the
SIMD kernels, which run at a much higher rate, the crossover
`FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF` is between 256 and 800,
depending on the processor. The byte kernels of `nmod_mat_mul_u8` recurse
from dimension 256 ($N \le 16$) or 2048 ($17 \le N \le 255$), and the
bit-packed kernels for $N = 2$ and $N = 3$ from 8192 and 16384, since
their leaves are much cheaper relative to an addition.

Numerical stability, the main concern with Strassen's algorithm in
floating point, plays no role over $\Z/N\Z$, where all operations are
exact.
