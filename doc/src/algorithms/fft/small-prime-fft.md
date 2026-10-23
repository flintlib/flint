---
title: FFTs modulo word-size primes in floating point
summary: The number-theoretic transforms of `fft_small` modulo primes $p < 2^{50}$, on integers stored in doubles: the order of the roots, radix-4 butterflies with lazy reduction, blocking and truncation, and proofs that every value stays in range.
complexity: A transform of length $2^L$ takes $\tfrac12 L 2^L$ modular multiplications by roots of unity, $\tfrac18 L 2^L$ reductions and $L 2^L$ additions, all on vectors of 8 doubles; a truncated transform costs roughly in proportion to the number of points computed.
order: 40
---

The module [`fft_small`](fft_small) computes number-theoretic transforms,
that is, discrete Fourier transforms over $\Z/p\Z$ [@Pol1971], for primes
$p < 2^{50}$ with $2^L \mid p - 1$. Residues are integers stored in
`double` variables and the modular arithmetic is done in floating point
with fused multiply-adds, four lanes per instruction with AVX2 and two
with NEON, the code always working on vectors of eight. This approach to
SIMD modular arithmetic is that of van der Hoeven, Lecerf and Quintin
[@HLQ2016]. The implementation, by Daniel Schultz, combines radix-4
butterflies with lazy reduction in the spirit of Harvey [@Har2014], van
der Hoeven's truncated transforms [@vdH2004] and a recursive matrix
decomposition of the transform [@Bai1990]. How several primes are combined is described in
[](algorithms/fft/multimodular), and how polynomial multiplication uses
the transforms in [](algorithms/polynomials/multiplication#fft-small).

Throughout, $p$ is the prime and $u = 2^{-53}$. Values are integers of
either sign and are never reduced to $[0, p)$ inside a transform; a
*bound* on a value means a bound on its absolute value.

## Primes and roots of unity {#roots}

A prime is usable if it passes the check of [](#arithmetic), which
implies $p < 2^{50}$, if $2^4 \mid p - 1$, and if $2^L \mid p - 1$ for
the largest transform length $2^L$ used. If $2^v \| p - 1$ and $a$ is a
quadratic nonresidue, $g = a^{(p-1)/2^v}$ has order exactly $2^v$,
because $g^{2^{v-1}} = a^{(p-1)/2} = -1$. All roots of unity are powers
of this $g$.

The transforms evaluate at the roots of unity in *bit-reversed order*.
For $j < 2^k$ let $\operatorname{rev}_k(j)$ be $j$ with its $k$ low bits
reversed, and let $\zeta_{2^k} = g^{2^{v-k}}$, so that
$\zeta_{2^{k+1}}^2 = \zeta_{2^k}$. Define

$$
w_j = \zeta_{2^k}^{\operatorname{rev}_k(j)} \qquad (j < 2^k).
$$

This does not depend on $k$, since $\operatorname{rev}_{k+1}(j) =
2 \operatorname{rev}_k(j)$ for $j < 2^k$. The sequence begins
$1, -1, \zeta_4, -\zeta_4, \zeta_8, \zeta_8^5, \zeta_8^3, \zeta_8^7,
\ldots$; `sd_fft_ctx_w` returns its terms.

::: lemma {#lem-w}
For every $j \ge 0$: $w_{2j}^2 = w_j$ and $w_{2j+1} = -w_{2j}$.
Consequently, for every $h \ge 0$,

$$
x^{2^h} - w_j = \prod_{i < 2^h} \left(x - w_{j 2^h + i}\right).
$$
:::

::: proof
Let $j < 2^k$. The bits of $2j$ are those of $j$ shifted up, so
$\operatorname{rev}_{k+1}(2j) = \operatorname{rev}_k(j)$ and
$w_{2j}^2 = \zeta_{2^{k+1}}^{2 \operatorname{rev}_k(j)} = w_j$. The low
bit of $2j + 1$ becomes the top bit of the reversal, so
$\operatorname{rev}_{k+1}(2j+1) = 2^k + \operatorname{rev}_k(j)$ and
$w_{2j+1} = \zeta_{2^{k+1}}^{2^k} w_{2j} = -w_{2j}$. For the product,
use induction on $h$: by the first part
$x^{2^h} - w_j = (x^{2^{h-1}} - w_{2j})(x^{2^{h-1}} - w_{2j+1})$, and the
indices $2j \cdot 2^{h-1} + i$ and $(2j + 1) 2^{h-1} + i$ for
$i < 2^{h-1}$ are exactly the $j 2^h + i$ for $i < 2^h$.
:::

Since $w_{2j+1} = -w_{2j}$, only the even-indexed roots are stored, as
integers in $[-p/2, p/2]$. Table $k \ge 2$ holds $w_{2j}$ for
$2^{k-1} \le j < 2^k$, and it is obtained by multiplying all earlier
stored roots by $\zeta_{2^{k+1}}$, since $w_{2j} = \zeta_{2^{k+1}}
w_{2j - 2^k}$ in that range. Tables are added on demand under a mutex and
never move, so the transforms read them without locking. The inverse
transform needs $-w_{2j}^{-1}$; since $w_{2j}^{-1}$ is again a root of
the same order, it is found in the same table at the mirrored index
(`sd_fft_ctx_w2inv`).

## The forward transform {#forward}

By [](#lem-w), the values $f(w_i)$ for $i < 2^L$ can be computed by a
tree of remainders [@Ber2008]. At level $h$, block $j$ (for
$j < 2^{L-h}$) holds the $2^h$ coefficients of $f \bmod (x^{2^h} - w_j)$;
at level $L$ there is the single block $f \bmod (x^{2^L} - 1) = f$, and
at level 0 block $i$ holds $f \bmod (x - w_i) = f(w_i)$.

::: algorithm Radix-2 step {#alg-radix2}
**Input:** block $j$ at level $h$, with lower half $a$ and upper half $b$
(the residue is $\sum_i a_i x^i + x^{2^{h-1}} \sum_i b_i x^i$), and
$s = w_{2j}$.

**Output:** blocks $2j$ and $2j + 1$ at level $h - 1$.

1. For each $i < 2^{h-1}$: $c_i \gets a_i + s b_i$, $d_i \gets a_i - s b_i$.
2. Block $2j$ is $c$, block $2j + 1$ is $d$.
:::

This is correct because $x^{2^{h-1}} \equiv s$ modulo
$x^{2^{h-1}} - w_{2j}$ and $x^{2^{h-1}} \equiv -s$ modulo
$x^{2^{h-1}} - w_{2j+1}$. The multiplications by $s$ are the only modular
multiplications, and they vanish for $j = 0$, where $s = 1$. Two
consecutive levels are merged into *radix-4 butterflies*: with the
quarters $a_0, a_1, a_2, a_3$ of block $j$ and $s_1 = w_{2j}$,
$s_2 = w_{4j}$, $s_3 = w_{4j+2} = \pm \zeta_4 s_2$,

$$
\begin{aligned}
b_0 &= (a_0 + s_1 a_2) + s_2 (a_1 + s_1 a_3), &
b_1 &= (a_0 + s_1 a_2) - s_2 (a_1 + s_1 a_3), \\
b_2 &= (a_0 - s_1 a_2) + s_3 (a_1 - s_1 a_3), &
b_3 &= (a_0 - s_1 a_2) - s_3 (a_1 - s_1 a_3),
\end{aligned}
$$

which costs four multiplications for four outputs and keeps the data in
registers for two levels. A transform has $\lfloor L/2 \rfloor$ such
levels, plus a radix-2 level or an 8-point kernel when $L$ is odd. The
output order is the order of the blocks, so position $i$ holds
$f(w_i)$, except that the innermost kernel saves its final $4 \times 4$
transposition: for $L \ge 4$ the value $f(w_i)$ is at position
`sd_fft_ctx_trunc_index(L, i)`, obtained by exchanging the two 2-bit
fields of the low four bits of $i$. Pointwise products do not depend on
the order, and the inverse transform undoes it.

### Blocking {#blocking}

The data is divided into blocks of 256 doubles (2 KB). Up to 256 points,
the transform is done by kernels that keep the data in registers: radix-4
levels on vectors of eight lanes down to 16 points, then $4 \times 4$
transposes and two levels on 4-lane vectors. A longer transform of $2^k$
blocks treats the blocks as the entries of a $2^{k_1} \times 2^{k_2}$
matrix, $k_1 = \lfloor k/2 \rfloor$, stored by rows. The first $k_1$
levels are transforms of length $2^{k_1}$ along the columns, each acting
on whole blocks at stride $2^{k_2}$; the remaining levels are independent
transforms of the rows, each a contiguous transform of $2^{k_2}$ blocks
with the roots of the row's position $j 2^{k_1} + b$ in the tree. Both are
done recursively in the same way. This is Bailey's four-step method
[@Bai1990] in the remainder-tree formulation, which needs neither
twiddle multiplications between the steps nor a transposition (compare
the formulation with twiddle factors in
[](algorithms/fft/matrix-fourier)), and the
recursive halving makes it cache-oblivious in the sense of [@FLPR1999]:
every subproblem eventually fits in some level of the cache. In the
outer levels each butterfly is a loop over the 256 lanes of a block,
which gives long vectorised loops with one set of roots.

### Truncation {#truncation}

For a product of length $t < 2^L$, only $t$ values are needed and the
inputs have at most $t$ nonzero coefficients. `sd_fft_trunc` takes an
input truncation (coefficients from there on are zero) and an output
truncation (only positions below it are needed), in whole blocks for
$L > 8$. In the matrix decomposition, a column transform is skipped if
all its inputs are zero and computes only as many outputs as there are
needed rows, and of the row transforms only the needed ones are
computed, the last one truncated. At the radix-4 level, butterflies with
known zero inputs or unneeded outputs are replaced by special versions,
one for each of the 12 combinations of 2 to 4 inputs and 1 to 4 outputs.
Since a product of length $t \in (2^{L-1}, 2^L]$ needs a transform of
length $2^L$, truncation avoids the factor of up to two that an
untruncated transform would cost [@vdH2004].

## The inverse transform {#inverse}

The inverse runs the tree from the leaves to the root. If block $2j$
holds $c = A + sB$ and block $2j+1$ holds $d = A - sB$ (all operations
coefficient-wise), with $s = w_{2j}$, then $c + d = 2A$ and
$(c - d) s^{-1} = 2B$. The code computes $(c + d, (d - c) W)$ with
$W = -s^{-1}$, in radix-4 form, like the butterflies of Gentleman and
Sande [@GS1966]. The factors 2 are not removed, so the result is $2^L$
times the input polynomial; the factor $2^{-L}$ is applied together with
the pointwise product.

The truncated inverse [@vdH2004] starts from the values at the first $t$
positions of the transform of a polynomial of length at most $t$. In a
block where some values and some coefficients are known, the unknowns
are solved for from $c = A + sB$ and $d = A - sB$ instead of inverting
the butterfly: if $c$ and $B$ are known, then $A = c - sB$ and the value
$d = c - 2sB$ follows, which the next level up may need. The recursion of
`sd_ifft_trunc` follows the matrix decomposition of the forward transform,
with three parameters per subproblem (the number of known entries, the
number of wanted coefficients, and whether one further value is wanted),
and special radix-2 and radix-4 butterflies for the combinations that
occur. Truncations are in whole blocks for $L > 8$, and for $L \le 8$ the
full transform is computed. See also [](algorithms/fft/truncated-fft).

## Transposed transforms {#transposed}

`sd_fft_trunc_tr` and `sd_ifft_trunc_tr` compute the transposes of the
linear maps $d \mapsto$ `sd_fft_trunc`$(d)$ and $d \mapsto$
`sd_ifft_trunc`$(d)$, with the same parameters, so that
$\langle F(a), v \rangle = \langle a, F^T(v) \rangle$ for every input
pair. They are used for the transposed products of
[](algorithms/polynomials/middle-product#transposed-fft), which compute
middle products at the cost of a product of the output length. Their
derivation is mechanical. The transforms are straight-line programs of
butterflies, block copies and block zeroings, and the transpose of such a
program performs the transposed steps in the reverse order
([](algorithms/polynomials/middle-product#thm-tellegen)):

- The transpose of the forward butterfly $(a, b) \mapsto (a + wb,
  a - wb)$ is $(a, b) \mapsto (a + b, w(a - b))$: a butterfly of the
  inverse (Gentleman–Sande) shape with the forward twiddle $w$. So
  `sd_fft_trunc_tr` uses the twiddle tables of `sd_fft_trunc` and the
  reduction pattern of `sd_ifft_trunc`, and conversely
  `sd_ifft_trunc_tr` uses forward butterflies with the inverse twiddles.
  In the source the four transforms share one set of butterfly macros
  (`src/fft_small/butterflies.h`), which take the modular multiplication
  as an argument: the inverse-shape butterfly computes
  $(a, b) \mapsto (a + b, W(b - a))$, and the transposed transforms
  pass a negated multiplication, so that $W = -w$ without negating the
  twiddle tables. Every intermediate value then has the magnitude it has
  in the transform whose butterflies are shared, and the bound analyses
  below apply unchanged.
- The recursion is traversed in reverse: in the matrix decomposition the
  row transforms come before the column transforms, and within the
  basecases the levels run from last to first. The $4 \times 4$
  transposition that the forward kernels omit at the end (see
  `sd_fft_ctx_trunc_index`) becomes an omitted transposition at the start,
  so the transposed forward transform reads its input in the order in
  which the forward transform writes its output.
- Truncation transposes into truncation: not reading input $i$ deletes a
  column of the matrix and not writing output $o$ deletes a row, so the
  transposed truncated forward butterflies are the same macros with loads
  guarded by `otrunc` and stores by `itrunc`. The truncated inverse
  butterflies of [](#inverse) are separately derived maps (thirteen
  radix-4 and six radix-2 variants), and each is transposed separately.

Each transposed transform therefore performs exactly the operations of
the transform it transposes, and costs the same. `sd_fft_trunc_tr` reads
the evaluation side and accepts the range $(-2p, 2p)$ of `sd_ifft_trunc`
inputs, which pointwise products produce, because its butterflies have
the inverse shape; it also returns values in $(-2p, 2p)$, which was
checked by interval evaluation of its truncated butterflies and by
propagating bounds through whole transforms up to depth 16.
`sd_ifft_trunc_tr` has the forward shape and maps $(-3p, 3p)$ into
itself, as `sd_fft_trunc` does ([](#thm-forward)). The tests check the
adjoint identity on random vectors over a range of depths and
truncations.

## Pointwise products and normalisation {#pointwise}

After the forward transforms of $a$ and $b$, the pointwise product and
the inverse transform give $2^L$ times the cyclic convolution of $a$ and
$b$ of length $2^L$, which is their product if its length is at most the
truncation. `sd_fft_ctx_point_mul` computes $m x_i y_i$ with a constant
$m$ in the same pass, and the multiplication drivers take
$m = (c 2^L)^{-1} \bmod p$, where $c$ is the reduced Chinese remainder
cofactor of $p$, so that neither the normalisation nor the inverse in the
Chinese remainder formula costs anything extra.

## Arithmetic and its ranges {#arithmetic}

Each step is one of the following operations of
[](algorithms/simd/fp-modular-arithmetic), applied to integers stored
exactly in doubles:

- `reduce_to_pm1n` maps an integer $a$ with $|a| < 2^{53}$ to an integer
  $r \equiv a \pmod p$ with $|r| < p/2 + 2$
  (by [the reduction theorem](algorithms/simd/fp-modular-arithmetic#thm-reduce));
- `mulmod_fast(a, w)` and `mulmod(a, w)` return $r \equiv aw \pmod p$;
- additions and subtractions, which are exact for results below $2^{53}$,
  in particular for all sums below $8p$.

For the products FLINT checks each prime with
`fft_small_mulmod_satisfies_bounds`, from `sd_fft_ctx_init_prime`, for
the two properties

$$
\text{(P1)} \quad |aw| < 2p^2 \implies |r| < p, \qquad
\text{(P2)} \quad |aw| < 4p^2 \implies |r| < \tfrac32 p,
$$

(P1) for both functions and (P2) for `mulmod`. The check evaluates the
error bound of [the multiplication theorem](algorithms/simd/fp-modular-arithmetic#thm-mulmod)
with the actual error of $\mathrm{fl}(1/p)$ instead of its maximum; every
modulus below $2^{49}$ passes, and about one 50-bit prime in 900 fails.
The roots of unity satisfy $|w| \le (p - 1)/2$.

::: lemma Forward butterflies {#lem-forward}
Let $p \ge 16$ satisfy (P1). If the inputs of a forward butterfly (radix
2, radix 4, the 8-point kernels, or a truncated version) are less than
$3p$, then every product passed to `mulmod_fast` is less than $2p^2$,
every sum is less than $6p$, and the outputs are less than
$\tfrac52 p + 2$.
:::

::: proof
Consider the radix-4 butterfly for $j \ne 0$. The code reduces $a_0$ to
$|a_0'| < p/2 + 2$, and multiplies $a_2$ and $a_3$ by $s_1$, with
products below $3p \cdot p/2 < 2p^2$, so the results $a_2', a_3'$ are
less than $p$ by (P1). Then $a_0' \pm a_2'$ is less than
$\tfrac32 p + 2$, and $a_1 \pm a_3'$, with $a_1$ unreduced, less than
$4p$. The latter are multiplied by $s_2$ or $s_3$, with products below
$4p \cdot p/2 = 2p^2$, giving results less than $p$. The outputs are
sums of a value below $\tfrac32 p + 2$ and one below $p$. For $j = 0$ the
butterfly reduces $a_0$, $a_2$, $a_3$ and the sum $a_1 + a_3$ instead of
multiplying by $s_1 = s_2 = 1$, and multiplies $a_1 - a_3$ (below
$\tfrac72 p + 2$) by $s_3 = \zeta_4$, with product below $2p^2$; its
outputs are below $2p + 4$. The radix-2 butterfly gives outputs below
$\tfrac32 p + 2$. The 8-point kernel for $j \ne 0$ performs three radix-2
levels, with a reduction or a multiplication applied to every value
before the last additions, so its outputs are below $\tfrac32 p + 2$; the
products it forms have one factor below $4p$ and one at most $p/2$. The 8-point kernel for
$j = 0$ reduces all eight sums of its first level, and its outputs are
sums of at most four values below $p/2 + 2$, or of one or two such values
and a product, all below $\tfrac52 p + 2$. The truncated versions perform
the same operations with some inputs replaced by zero, or copy an input
to several outputs, which does not increase any bound.
:::

::: theorem Forward transform {#thm-forward}
Let $p \ge 16$ satisfy (P1). If the inputs of `sd_fft_trunc` are integers
of absolute value less than $A \le 3p$, the computation is exact and the
outputs are integers less than $\max(A, \tfrac52 p + 2) \le 3p$ in
absolute value, congruent to the values of the transform.
:::

::: proof
Every value written is either an input, a copy of a value written before
(in a truncated block with a single nonzero input), or an output of a
butterfly whose inputs were written before. By induction and
[](#lem-forward), all values are less than $\max(A, \tfrac52 p + 2)$;
all floating-point operations act on integers below $2^{53}$ and the
modular operations satisfy their hypotheses, so every result is the exact
integer, congruent to the corresponding coefficient of the remainder
tree.
:::

The conversions of the multiplication drivers produce inputs below
$p + 2^{32}$, so the outputs of their forward transforms are below
$\tfrac52 p + 2$.

::: proposition Pointwise products {#prop-pointwise}
Let $0 \le m \le p$. For $|x|, |y| < 3p$, `sd_fft_ctx_point_mul` returns
a value $r \equiv m x y \pmod p$ with $|r| < \tfrac32 p$. For
$|x| < \tfrac52 p + 2$, `sd_fft_ctx_point_sqr` returns
$r \equiv m x^2 \pmod p$ with $|r| < p$. Both are exact.
:::

::: proof
The constant is first brought to $|m'| \le p/2$. For the product,
$|x m'| < \tfrac32 p^2 < 2p^2$ gives a result $t$ with $|t| < p$ by (P1),
and $|t y| < 3p^2 < 4p^2$ gives $|r| < \tfrac32 p$ by (P2); pointwise
products use `mulmod`. The square $x^2$ is formed first:
$x^2 < 6.26 p^2 < 2^{103}$, which `mulmod` reduces exactly to a value below
$p/2 + 3.0001 \cdot 2^{-53} \cdot 6.26 p^2 < 2.85 p$
(by [the multiplication theorem](algorithms/simd/fp-modular-arithmetic#thm-mulmod)). Its product with
$m'$ is below $1.43 p^2 < 2p^2$, and (P1) gives $|r| < p$.
:::

::: lemma Inverse butterflies {#lem-inverse}
Let $p$ satisfy (P1) and (P2). If the inputs of an untruncated inverse
butterfly (radix 2, radix 4, or the 2-, 4- and 8-point kernels) are less
than $2p$, then the products passed to `mulmod_fast` are less than
$2p^2$, those passed to `mulmod` less than $4p^2$, all sums less than
$8p$, and the outputs less than $2p$.
:::

::: proof
For the radix-4 butterfly with $j \ne 0$, the sums $y_0 = x_0 + x_1$ and
$y_1 = x_2 + x_3$ are below $4p$, and the differences $x_0 - x_1$ and
$x_3 - x_2$, also below $4p$, are multiplied by roots of absolute value
at most $p/2$, with products below $2p^2$ and results below $p$. The
outputs are: $y_0 + y_1$ (below $8p$) reduced to below $p/2 + 2$; the
difference of the two products, below $2p$; $y_1 - y_0$ (below $8p$)
multiplied by a root with `mulmod`, product below $4p^2$ and result below
$\tfrac32 p$ by (P2); and the sum of the two products (below $2p$)
multiplied by a root, result below $p$. For $j = 0$, three of the sums
are reduced and the fourth multiplied by a fourth root of unity (product
below $2p^2$), and the outputs, sums of two of these values, are below
$\tfrac32 p + 2 < 2p$. The radix-2 butterflies return a reduced sum and a
product below $p$, and the kernels are compositions of these butterflies.
:::

So an untruncated inverse transform maps $(-2p, 2p)$ into itself, and by
[](#prop-pointwise) the pointwise products are valid inputs.

::: remark Truncated inverse transform {#rem-truncated-inverse}
For the truncated inverse butterflies the analogue of [](#lem-inverse)
does not follow from (P1), (P2) and the bound $p$ on reductions: five of
them can then produce values up to $5p$ from inputs below $2p$, and
through the recursion such values would reach `mulmod_fast` with products
up to $3p^2$. `sd_ifft_trunc` is nevertheless correct, because the actual
ranges of the operations are smaller: `reduce_to_pm1n` returns less than
$p/2 + 2$ and the products less than $p/2 + 3u|aw|$. This was verified
mechanically by running `sd_ifft_trunc` on bounds instead of values: each
entry holds a bound on the absolute value it stands for, each operation
checks its hypothesis and returns the smaller of its guaranteed and its
actual range, all roots of unity other than $\pm 1$ are replaced by the
bound $(p - 1)/2$, the error of $\mathrm{fl}(1/p)$ by its maximum, and $p$
by a value close to $2^{50}$, which maximises the rounding errors relative
to $p$. For inputs below $2p - 1$, all truncations, depths up to 25 and
both rounding methods, every product given to `mulmod_fast` is below
$2p^2$, every sum below $8p$, the values passed between butterflies are
below $3.5 p$ and the outputs at most $2.5 p$. Outputs above $2p$ do
occur, so the stated output range of `sd_ifft_trunc` is $(-3p, 3p)$, the
input range of `sd_fft_trunc`; FLINT's callers reduce the outputs or
multiply them by a reduced constant first. The derivation is recorded in
`fft_small/mulmod_satisfies_bounds.c`; five more reductions would make the
truncated inverse provably map $(-2p, 2p)$ into itself, at a cost of
about 1%.
:::

::: proposition Two-prime pointwise products {#prop-two-prime}
The two-prime integer multiplication (`_mpn_ctx_mpn_mul_window_np2`)
multiplies two forward outputs $x, y$ with `mulmod` without reducing
either first. For the two primes it uses, $p_0 = 63 \cdot 2^{44} + 1$ and
$p_1 = 75 \cdot 2^{43} + 1$, the results are exact and below $1.70 p_0$
and $1.26 p_1$ respectively, valid inputs for the inverse transform, with
either rounding method.
:::

::: proof
The inputs of the forward transforms are chunks of at most 45 bits, so
$|x|, |y| < \tfrac52 p + 2$ by [](#thm-forward), and $|xy| < 6.26 p^2 <
2^{103}$, outside the hypotheses of (P1) and (P2). Following the analysis
behind `fft_small_mulmod_satisfies_bounds`, the computed quotient $q$
satisfies

$$
\left| \frac{xy}{p} - q \right| \le \frac{|xy| \, t}{p}
  + \frac{|\ell|}{p} + \frac12 + 2^{-(B+1)},
$$

where $t = |p \, \mathrm{fl}(1/p) - 1|$, $|\ell| \le 2^{e - 54}$ is the
low part of the product for $|xy| < 2^e$, and $B$ is the number of
fractional bits kept when the approximate quotient is rounded to a
double. For $p_0$: $t = 0.248 u$, $e = 103$, $B = 0$, and the three error
terms are $0.191$, $0.508$ and $1$ (counting the last two together), so
$|xy - q p_0| < 1.70 p_0$. For $p_1$: $t = 0.169 u$, $e = 102$, $B = 1$,
and the terms are $0.077$, $0.427$ and $0.75$, which gives $1.26 p_1$.
The result is exact as in
[the multiplication theorem](algorithms/simd/fp-modular-arithmetic#thm-mulmod), since
$|xy| < 2^{103}$.
:::

The guaranteed ranges alone would allow $|xy|$ up to $9p^2$ and give no
bound below $2p$, so this step depends on the two specific primes. All
other drivers multiply by the normalisation constant, or reduce, before
forming the product, as in [](#prop-pointwise).

## Cost {#cost}

A radix-4 butterfly with $j \ne 0$ performs four multiplications by
roots, one reduction and eight additions for four points and two levels,
so a full transform of length $2^L$ takes $\tfrac12 L 2^L$
multiplications, $\tfrac18 L 2^L$ reductions and $L 2^L$ additions, on
8-lane vectors. A multiplication costs two floating-point multiplications,
two fused multiply-adds, a rounding and an addition. Truncated transforms
cost roughly in proportion to the number of points computed. The inverse
transform has the same cost. The data is transformed in place, and the
twiddle tables take half a double per point of the longest transform.
The lengths from which these transforms beat the other polynomial
multiplication algorithms are given in
[](algorithms/polynomials/multiplication#cutoffs).
