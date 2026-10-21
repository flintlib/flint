---
title: Basecase multiplication kernels
summary: Schoolbook multiplication of short limb arrays by a table of straight-line routines, one for each pair of small sizes (assembly with two independent carry chains on x86-64 and ARM64, C elsewhere), and Karatsuba multiplication built on them.
complexity: $mn$ limb multiplications for an $m \times n$ product; $O(n^{\log_2 3})$ with Karatsuba above the table.
order: 12
---

For operands of up to a few dozen limbs, integer multiplication is
dominated not by the number of limb multiplications but by everything
around them: loop control, function calls, branch mispredictions, loads
and stores, and above all the latency of carry propagation. GMP's
`mpn_mul_basecase` is a loop over variable-length row operations
(`mpn_mul_1`, then `mpn_addmul_1` for each further limb of the shorter
operand), each of which pays a start-up cost and ends its carry chain at
the row boundary. Ahlbäck and Johansson [@AJ2025] replace it by a *table
of straight-line routines*, one for each pair of small sizes
([](#performance)).

FLINT's `flint_mpn_mul`, `flint_mpn_mul_n` and `flint_mpn_sqr`
therefore look up a routine in such a table for small sizes, and switch
to Karatsuba recursion over these routines before using GMP's
Toom–Cook code and finally the FFT
([](algorithms/integers/multiplication)). The routines are written in
assembly on x86-64 processors with the BMI2 and ADX extensions and on
ARM64. Other builds have portable C routines for at most 7 limbs (and
for a one-limb operand times up to 14 limbs), and call GMP's
`mpn_mul_basecase` directly when the longer operand has at most 16
limbs ([](#hardcoded)). The selection is given in [](mpn_extras/mul),
and the high and low products built with the same techniques in
[](algorithms/integers/short-products).

## Schoolbook multiplication {#schoolbook}

Let $B = 2^{\bits}$ and let $x = \sum_{i<m} x_i B^i$,
$y = \sum_{j<n} y_j B^j$. The schoolbook product computes

$$
xy = \sum_{j=0}^{n-1} \Bigl(\sum_{i=0}^{m-1} x_i y_j B^i\Bigr) B^j,
$$

one row $x y_j$ at a time (*operand scanning*), adding each row into the
partial result shifted by $j$ limbs. Each limb product $x_i y_j$ is a
two-limb number $(h, \ell)$ with $x_i y_j = hB + \ell$.

::: lemma {#lem-mac}
For limbs $a, b, c, d < B$, the value $ab + c + d$ is less than $B^2$, so
a multiply-accumulate of a limb product and two limbs never overflows two
limbs.
:::

::: proof
$ab + c + d \le (B-1)^2 + 2(B - 1) = B^2 - 1$.
:::

[](#lem-mac) is what makes the row operation $r \gets r + x y_j$ possible
with a single carry limb per position: at position $i$ one adds $\ell_i$
(the low half of $x_i y_j$), the high half $h_{i-1}$ of the previous
product and the incoming carry to $r_i$. Done naively, this needs two
additions with carry per limb that depend on each other through the carry
flag, so that a row runs at the latency of a chain of $2m$ dependent
additions.

::: algorithm Row multiply-accumulate with two carry chains {#alg-row}
**Input:** $x_0, \ldots, x_{m-1}$, a limb $y$, a partial result
$r_0, \ldots, r_{m-1}$.

**Output:** $r + xy$ in $r_0, \ldots, r_m$.

1. For $i = 0, \ldots, m - 1$: compute $(h_i, \ell_i) = x_i y$ without
   affecting the flags (`mulx`).
2. Chain A (one flag): add $\ell_i$ to $r_i$ for
   $i = 0, \ldots, m - 1$, propagating its carry.
3. Chain B (the other flag): with $r_m = 0$, add $h_i$ to
   $r_{i+1}$ for $i = 0, \ldots, m - 1$, propagating its carry.
4. Add the final carry of chain A into $r_m$.
:::

On x86-64 processors with the BMI2 and ADX extensions (Intel Broadwell and
later, AMD Zen and later), `mulx` multiplies without touching the flags
and `adcx` and `adox` add with carry through two different flags, so the
two chains of [](#alg-row) are interleaved and execute in parallel, and
the multiplications for the next limbs are issued while the additions are
in flight. Since adding a limb to a limb produces a carry of at most one,
each chain carries one bit per position. In the schoolbook product the
partial result in the window of row $j$ is less than $B^m$, so the row
result $r + x y_j \le (B^m - 1) + (B^m - 1)(B - 1) < B^{m+1}$ fits in
$m + 1$ limbs: chain B ends without a carry, and step 4 cannot
overflow.

In GMP, the rows are separate calls to `mpn_addmul_1`, a loop unrolled
a fixed number of times $k$, and each call has its own set-up and its own
carry chain. Ahlbäck and Johansson analyse the Apple M1, which can
issue one limb multiplication per cycle: there an `mpn_addmul_1` loop
unrolled $k$ times needs about $(k+1)/k$ cycles per limb product, GMP's
reaching the bound of 5 cycles per 4 limbs set by the units that handle
the carry flag [@AJ2025]. Writing the whole product as one routine
removes the calls, the loop control and the row boundaries, and lets
the processor overlap the carry chains of consecutive rows; the target
is one cycle per limb product.

## Hard-coded sizes {#hardcoded}

For small $m$ and $n$, the whole product is written out as straight-line
code: there is no loop and no branch, the running row sum stays in
registers where possible (on x86-64 each row multiplies the limbs of the
longer operand by one limb of the shorter, held in the implicit operand
register of `mulx`, and retires its lowest limb to memory), and the
carry chains of consecutive rows are scheduled to overlap. Each size
pair is a separate function, and the dispatching functions look the
routine up in a table indexed by the sizes. The two architectures call
for different shapes [@AJ2025]:

- **x86-64 with ADX.** The two-address instruction format and the 15
  usable general-purpose registers make it necessary to unroll
  completely, in both sizes, to avoid loads, stores and moves. The
  assembly routines `flint_mpn_mul_m_n` for all $1 \le n \le m$ with
  $m \le 16$ and $n \le 8$ (100 routines) are generated by a script
  (`dev/gen_mul_basecase.jl`). Larger cases of the table are C functions
  built from these: products with $9 \le n \le 14$ and $m \le 16$ use the
  routine for $n = 8$ followed by $n - 8$ calls to `mpn_addmul_1`
  (*peeling*); $15 \times 15$ is split into $15 \times 8$ and
  $15 \times 7$, and $16 \times 15$ and $16 \times 16$ use Karatsuba.
  Squares have their own routines up to 14 limbs, and high and low
  products their own (see [](algorithms/integers/short-products)).
- **ARM64.** The three-address format and 29 usable registers allow a
  *semi-unrolled* kernel: `flint_mpn_mul_kn`, for each $k \le 15$, keeps
  a $k$-limb operand in registers (from $k = 9$ on, parts of it are
  reloaded) and loops over the limbs of the other operand, of any length.
  `flint_mpn_mul_2` multiplies by two limbs; squares have routines up to
  9 limbs.
- **Other machines.** Portable C routines for $m \le 7$ (and $n = 1$ up to
  $m = 14$), generated by `dev/gen_mul_basecase.py`. They use product
  scanning: each output limb is the sum of a diagonal of limb products,
  accumulated in a three-limb accumulator held in local variables that
  the compiler keeps in registers. Above, GMP's `mpn_mul_basecase` is
  called up to 16 limbs.

The table sizes balance the speed of each routine against the size of
the code: very large unrolled routines run out of registers and fill the
instruction cache, and larger tables make the indirect jump of the
dispatch harder to predict. Above the hard-coded sizes, an $m \times n$
product can be split into hard-coded blocks by *peeling* (one large block
plus rows), by a *balanced* split into equal blocks, or by Karatsuba;
which is fastest for a given pair of sizes depends on the machine and is
found by measurement [@AJ2025].

In a square $x^2 = \sum_i x_i^2 B^{2i} + 2 \sum_{i<j} x_i x_j B^{i+j}$,
the off-diagonal products are computed once, doubled by a shift, and the
diagonal squares added, which needs about $n(n+1)/2$ limb products
instead of $n^2$. On x86-64, the squares of 1 to 10 limbs are assembly
routines; for $11 \le n \le 14$, $x = x_1 B^7 + x_0$ is split after 7
limbs and $x^2 = x_1^2 B^{14} + 2 x_0 x_1 B^7 + x_0^2$ is formed from two
hard-coded squares and one hard-coded $7 \times (n-7)$ product, doubled
and added in one `mpn_addmul_1` pass by the limb 2.

## Karatsuba above the table {#karatsuba}

Karatsuba's method [@KO1962] splits $x = x_1 B^k + x_0$ and
$y = y_1 B^k + y_0$ and uses

$$
xy = x_1 y_1 B^{2k} + \bigl(x_1 y_1 + x_0 y_0 - (x_1 - x_0)(y_1 - y_0)\bigr) B^k
     + x_0 y_0,
$$

three half-size products instead of four. With $M(n)$ the cost of an
$n \times n$ product, $M(n) = 3M(n/2) + O(n)$ gives
$M(n) = O(n^{\log_2 3}) = O(n^{1.585})$. GMP applies Karatsuba recursively
down to its own schoolbook routine; `flint_mpn_mul_toom22`, a port of
GMP's `mpn_toom22_mul`, recurses into `flint_mpn_mul` instead, so that the
leaves of the recursion are the hard-coded kernels. On x86-64 with ADX
this is faster than GMP up to 230 limbs (483 on ARM64, measured on an
Apple M1), for nearly balanced operands ($5n \ge 4m$); beyond, GMP's
Toom–Cook multiplication takes over. Builds without FLINT's assembly do
not use `flint_mpn_mul_toom22`; they call GMP's `mpn_mul` up to the FFT
threshold.

## Dispatch {#dispatch}

`flint_mpn_mul`, `flint_mpn_mul_n` and `flint_mpn_sqr` are inline
functions: a size test and an indirect call through the table, or a call
to the general code above the table. When the sizes are compile-time
constants, the compiler removes the test and the table lookup and calls
the kernel directly. When sizes vary unpredictably,
the indirect call is a source of branch mispredictions; in the
measurements of [@AJ2025], a table-based multiplication is still faster
than GMP on random sizes, but by less (a factor of about 1.8 for random
sizes up to 8 limbs on Zen 3, against 1.04 up to 64 limbs).

## Performance {#performance}

Ahlbäck and Johansson [@AJ2025] measure the kernels on an AMD Zen 3
(Ryzen 7 PRO 5850U at a fixed 1.9 GHz) and an Apple M1:

- Against GMP's `mpn_mul_basecase` on Zen 3, the hard-coded kernels are
  up to 2.2 times as fast (for $3 \times 3$ limbs), about 1.25 times as
  fast for $12 \times 12$ and 1.4 times for $16 \times 16$; the fully
  unrolled $m \times 8$ kernels need about 1.36 cycles per limb product.
- On the M1, where GMP has no native basecase and falls back to rows
  of `mpn_addmul_1`, the semi-unrolled kernels stay close to one cycle
  per limb product, between 1.3 and 4.3 times as fast as GMP.
- With Karatsuba over the kernels, `flint_mpn_mul` is between 1.14 (at
  $96 \times 96$ limbs) and 2.5 times as fast as GMP's `mpn_mul` on Zen 3,
  and between 1.25 and 3 times on the M1.

## Use in FLINT {#use}

The kernels are reached through `flint_mpn_mul`, `flint_mpn_mul_n`,
`flint_mpn_sqr` and `FLINT_MPN_MUL_WITH_SPECIAL_CASES`, which every
multiplication of `fmpz` integers, `mpn_mod` residues and `arf` and
`mp_real` mantissas goes through. The same techniques give the
fixed-size dot products and polynomial products of
[packed coefficients](mpn_extras#packed) and the high and low products of
[](algorithms/integers/short-products). The `nmod` and `fmpz_mod`
modules have their own kernels for one and two limbs.

::: functions flint_mpn_mul flint_mpn_mul_n flint_mpn_sqr flint_mpn_mul_toom22 :::
