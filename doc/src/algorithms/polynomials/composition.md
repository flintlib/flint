---
title: Composition and modular composition
summary: Computing $f(g)$, $f(g) \bmod x^n$ and $f(g) \bmod h$ by Horner's rule, divide and conquer, Taylor shifts, the baby-step giant-step method of Brent and Kung, and the algorithm of Kinoshita and Li, with the choices made by `nmod_poly`.
complexity: $f(g)$: $O(M(\ell_1 \ell_2) \log \ell_1)$ by divide and conquer. $f(g) \bmod x^n$: $O(\sqrt{n}\, M(n) + n^{(\omega+1)/2})$ by Brent–Kung and $O(M(n) \log n)$ by Kinoshita–Li. $f(g) \bmod h$, $\deg h = n$: $O(n M(n))$ by Horner and $O(\sqrt{n}\, M(n) + n^{(\omega + 1)/2})$ by Brent–Kung.
order: 60
---

Composition $f(g)$ of polynomials appears in three forms: exact
composition of polynomials, composition of power series truncated to
$x^n$ (which requires $g(0) = 0$), and modular composition
$f(g) \bmod h$, the basic operation of polynomial factoring over finite
fields [@GatSho1992; @KalSho1998]. In all three, Horner's rule
$f(g) = f_0 + g(f_1 + g(f_2 + \cdots))$ uses $\ell_1 - 1$
multiplications by $g$, which is too slow for long $f$. The faster
methods split $f$ into blocks and reuse powers of $g$, or (for power
series) reorganise the computation entirely. Brent and Kung
[@BrentKung1978] introduced the baby-step giant-step method, which is
also the basis of evaluating polynomials at matrices by Paterson and
Stockmeyer [@PS1973]. Kinoshita and Li [@KL2024] recently found a
quasi-linear algorithm for power series composition.

Throughout, $f$ has length $\ell_1$, $g$ has length $\ell_2$, $M(n)$ is
the cost of multiplying two polynomials of length $n$, and $\omega$ is an
exponent of matrix multiplication ($\omega \le \log_2 7$ by Strassen's
method).

## Generic version {#generic}

All the composition algorithms described here are implemented once in
`gr_poly` for any commutative ring, and `nmod_poly`, `fmpz_poly`,
`fmpz_mod_poly`, `fq_poly`, `arb_poly` and `acb_poly` call them with
their contexts. The generic dispatchers choose by size only, as
described below for `nmod_poly` (over balls, the cutoff between
Brent–Kung and Kinoshita–Li depends on the precision). The Taylor
shift, a ring method, uses Horner's scheme up to length 20 and divide
and conquer above, with `fmpz` and the ball rings supplying their own
versions. Modular
composition with the Brent-Kung algorithm multiplies matrices with
`gr_mat_mul`, so it uses the fastest matrix multiplication of the ring
([](algorithms/linear-algebra/overview)). Series reversion by Lagrange
inversion divides by integers and fails in small characteristic, where
Newton iteration is used. The same algorithms over other rings are
summarised in [](algorithms/polynomials/overview#composition).

## Composition of polynomials {#compose}

The result $f(g)$ has length $(\ell_1 - 1)(\ell_2 - 1) + 1$. Horner's rule
multiplies polynomials of growing length $i(\ell_2 - 1) + 1$ by $g$, which
costs about as much as $\tfrac{1}{2}\ell_1^2$ products of length $\ell_2$. Divide and
conquer instead pairs up the coefficients:

$$
f(g) = \sum_{i} (f_{2i} + f_{2i+1} g)\, (g^2)^i ,
$$

which is a polynomial of half the length in $g^2$, with polynomial
"coefficients" $f_{2i} + f_{2i+1} g$. Repeating with $g^2, g^4, \ldots$
gives:

::: algorithm Divide-and-conquer composition {#alg-dc}
**Input:** $f$, $g$.

**Output:** $f(g)$.

1. $h_i \gets f_{2i} + f_{2i+1} g$ for $0 \le i < \lceil \ell_1 / 2 \rceil$;
   $G \gets g^2$.
2. While there are more than two $h_i$: replace them by
   $h_{2i} + h_{2i+1} G$, and $G$ by $G^2$.
3. Return $h_0 + h_1 G$.
:::

::: proposition {#prop-dc}
[](#alg-dc) uses $O(M(\ell_1 \ell_2) \log \ell_1)$ operations.
:::

::: proof
In round $j$ there are about $\ell_1/2^j$ products of $G = g^{2^j}$, of
degree $2^j (\ell_2 - 1)$, with polynomials of about the same degree,
costing $O((\ell_1/2^j) M(2^j \ell_2)) = O(M(\ell_1 \ell_2))$, plus one
squaring of the same size. There are $O(\log \ell_1)$ rounds.
:::

### Taylor shifts {#taylor-shift}

The special case $g = x + c$, the Taylor shift $f(x + c)$, has a
quasi-linear algorithm [@ASU1975]; von zur Gathen and Gerhard [@GG1997]
compare this and other methods.

::: lemma {#lem-taylor}
If $1, 2, \ldots, \ell_1 - 1$ are invertible, then
$$
f(x + c) = \sum_{k} \frac{x^k}{k!} \sum_{j \ge k} (j!\, f_j)
\frac{c^{j-k}}{(j-k)!},
$$
so the coefficients of $k!\,f(x + c)$ form a correlation of the vectors
$(j! f_j)$ and $(c^i / i!)$, that is one product of length $\ell_1$.
:::

::: proof
Expand $(x + c)^j = \sum_k \binom{j}{k} x^k c^{j-k}$ and write
$\binom{j}{k} = j!/(k!(j-k)!)$.
:::

More generally, $g = a x^d + c$ is a Taylor shift by $c$, a scaling of the
coefficients by powers of $a$ and the substitution $x \to x^d$.

**In FLINT.** `nmod_poly_compose` uses Horner's rule
(`nmod_poly_compose_horner`) for $\ell_1 \le 7$. Otherwise it calls
`_gr_poly_compose`, which handles $g = a x^d + c$ by a Taylor shift,
scaling and inflation, uses Horner's rule for $\ell_1 \le 7$, and
[](#alg-dc) otherwise. `nmod_poly_taylor_shift` uses Horner's rule
($O(\ell_1^2)$ additions and multiplications) for lengths below 100, for
lengths exceeding the modulus and for $c = \pm 1$ below length 1000, and
[](#lem-taylor) (`nmod_poly_taylor_shift_convolution`) otherwise. See
also [](nmod_poly/compose).

## Composition of power series {#series}

Let $g(0) = 0$ and $n$ the target precision. Horner's rule and divide and
conquer work as before, truncating every product to length $n$. Both are
fast when the exact composition is not much longer than $n$, but in
general they need $\Theta(\ell_1)$ products of length $n$.

### Brent–Kung {#brent-kung-series}

Let $m \approx \sqrt{n}$ and split $f$ into blocks of $m$ coefficients:

$$
f(g) = \sum_{i} F_i(g)\, (g^m)^i, \qquad F_i = \sum_{j < m} f_{im + j} x^j .
$$

The *baby steps* are $g^0, g^1, \ldots, g^{m-1} \bmod x^n$. Then all
$F_i(g)$ together are the rows of a matrix product: the matrix of
coefficients of the $F_i$ (about $\ell_1/m$ rows and $m$ columns) times
the matrix whose rows are the baby steps ($m$ rows, $n$ columns). The
*giant steps* combine the $F_i(g)$ by Horner's rule in $g^m$.

::: theorem Brent–Kung [@BrentKung1978] {#thm-bk-series}
For $\ell_1 \le n$, $f(g) \bmod x^n$ can be computed with
$O(\sqrt{n}\, M(n) + n^{(\omega + 1)/2})$ operations.
:::

::: proof
The baby steps and the giant steps are $O(\sqrt{n})$ products modulo
$x^n$. The matrix product of size $\sqrt{n} \times \sqrt{n}$ by
$\sqrt{n} \times n$ splits into $\sqrt{n}$ products of square matrices of
size $\sqrt{n}$, which cost $O(\sqrt{n} \cdot n^{\omega/2})$.
:::

### Kinoshita–Li {#kinoshita-li}

The algorithm of Kinoshita and Li [@KL2024] rests on the identity, valid
for $\ell_1 \le n$,

$$
f(g) \bmod x^n = [y^{n-1}] \frac{P(y)}{1 - y\, g(x)} \bmod x^n, \qquad
P(y) = \sum_{i < n} f_i y^{n-1-i},
$$

which follows by expanding $1/(1 - yg) = \sum_j y^j g^j$. Writing
$Q(x, y) = 1 - y g(x)$, the bivariate rational function is simplified by
Graeffe steps: $Q(x, y) Q(-x, y) = V(x^2, y)$, so

$$
\frac{P}{Q(x, y)} = \frac{P \cdot Q(-x, y)}{V(x^2, y)} ,
$$

which halves the $x$-degree that matters while the $y$-degree of the
denominator doubles. After $O(\log n)$ steps the denominator no longer
depends on $x$; the problem is then a univariate power series division,
and the result is lifted back through the stored $Q$'s, keeping at each
level only the slice of $y$-coefficients that contributes to
$[y^{n-1}]$. Each level costs $O(M(n))$ with bivariate products done by
Kronecker substitution, for $O(M(n) \log n)$ in total.

**In FLINT.** `nmod_poly_compose_series` calls `_gr_poly_compose_series`,
which uses direct formulas for $n \le 3$ and for monomials $g = c x^d$,
Horner's rule if $\ell_1 \le 7$ or $n \le 7$, divide and conquer if
$\ell_1^2 < n$ or $(\ell_1 - 1)(\ell_2 - 1) + 1 < 4n$ (when the exact
composition is not much longer than $n$), Brent–Kung for $n < 500$ and
Kinoshita–Li (`_gr_poly_compose_series_kinoshita_li`) for $n \ge 500$.

**Reversion.** The compositional inverse of $f$ (with $f(0) = 0$ and
$f'(0)$ a unit), `nmod_poly_revert_series`, uses Johansson's fast
Lagrange inversion [@Joh2015b], a baby-step giant-step method, for
$n < 2500$, and Newton iteration on $f(g) = x$, which needs one
composition and one series division per step, for larger $n$ or when
Lagrange inversion fails because some $k < n$ is not invertible.

## Modular composition {#modular}

Let $h$ have degree $n$ with invertible leading coefficient and let
$f$, $g$ have degree less than $n$. Horner's rule computes
$f(g) \bmod h$ with $\ell_1 - 1$ multiplications modulo $h$, each costing
$O(M(n))$ with Newton division and a precomputed inverse of
$\operatorname{rev}(h)$ ([](algorithms/polynomials/newton#divrem)): $O(n M(n))$
in total. The Brent–Kung method works modulo $h$ just as modulo $x^n$:

::: algorithm Brent–Kung modular composition {#alg-bk}
**Input:** $f$ of length $\ell_1 \le n$, $g$ of degree less than $n$, $h$
of degree $n$.

**Output:** $f(g) \bmod h$.

1. $m \gets \lfloor \sqrt{n} \rfloor + 1$.
2. Baby steps: $A_i \gets g^i \bmod h$ for $0 \le i < m$ (an
   $m \times n$ matrix $A$ of coefficients).
3. Let $B$ be the $m \times m$ matrix whose row $i$ holds coefficients
   $im, \ldots, im + m - 1$ of $f$. Compute $C \gets BA$; row $i$ of $C$
   is $F_i(g) \bmod h$.
4. Giant steps: $H \gets g^m \bmod h$; return
   $\sum_i C_i H^i \bmod h$, evaluated by Horner's rule in $H$.
:::

Correctness follows from $f(g) = \sum_i F_i(g) (g^m)^i$ and the fact that
reduction modulo $h$ commutes with sums and products. Since
$m^2 > n \ge \ell_1$, the matrix $B$ holds all coefficients of $f$.

::: theorem {#thm-bk}
[](#alg-bk) uses $O(\sqrt{n}\, M(n) + n^{(\omega+1)/2})$ operations.
:::

::: proof
Steps 2 and 4 are $O(m)$ multiplications modulo $h$. Step 3 multiplies an
$m \times m$ by an $m \times n$ matrix, that is $n/m$ products of
$m \times m$ matrices, costing $O(n m^{\omega - 1}) = O(n^{(\omega+1)/2})$.
:::

With classical matrix multiplication the matrix step costs $n^2$
operations, as much as schoolbook arithmetic; the method pays off
because matrix multiplication over $\Z/p\Z$ is very efficient in
practice (Strassen, BLAS and delayed reduction, see
[](algorithms/linear-algebra/matrix-multiplication)), and because the
matrix $A$ depends only on $g$ and $h$, so it can be reused for many
$f$. Nüsken and Ziegler [@NZ2004] improve the exponent using rectangular
matrix multiplication, and Kedlaya and Umans [@KU2011] give an algorithm
of nearly linear complexity over finite fields; neither is implemented
in FLINT.

FLINT improves step 4 with rectangular splitting when the matrix has more
than 10 rows: it precomputes $H^0, \ldots, H^t$ for $t \approx \sqrt{m}$
(at most 30), accumulates each block $\sum_{j < t} C_{it + j} H^j$
without reduction (each term is a product of two polynomials of degree
less than $n$), reduces once per block, and combines the blocks by
Horner's rule in $H^t$. This replaces $m$ modular multiplications by
about $m$ plain products and $m/t + t$ reductions.

**In FLINT.**

- `nmod_poly_compose_mod` reduces $g$ modulo $h$ and uses Horner's rule
  (`nmod_poly_compose_mod_horner`) if $\deg h < 7$ or $\ell_1 \ge n + 1$,
  and Brent–Kung (`nmod_poly_compose_mod_brent_kung`, which requires
  $\ell_1 \le n$) otherwise. The latter calls `_gr_poly_compose_mod_brent_kung`,
  which computes the baby steps as $g^i = g^{\lceil i/2 \rceil} g^{\lfloor i/2 \rfloor}$
  and multiplies the matrices with `gr_mat_mul`.
- `nmod_poly_compose_mod_brent_kung_preinv` takes the inverse of the
  reversal of $h$ as an extra argument; it computes the baby steps with
  `_nmod_poly_powers_mod_preinv_naive`, the matrix product with
  `nmod_mat_mul`, $H = g^{\lfloor m/2 \rfloor} g^{m - \lfloor m/2 \rfloor}$
  and the giant steps with the rectangular splitting described above
  (`_nmod_poly_mod_matrix_rows_evaluate`).
- `nmod_poly_precompute_matrix` computes the matrix $A$ of baby steps
  once, and `nmod_poly_compose_mod_brent_kung_precomp_preinv` reuses it
  for several $f$. The distinct-degree and equal-degree factorisation and
  irreducibility tests of `gr_poly` use the same technique to compose many
  polynomials with the same $g$ modulo the same $h$ [@KalSho1998].
- `nmod_poly_powers_mod_bsgs` computes $f^0, \ldots, f^{k-1} \bmod g$ by a
  baby-step giant-step decomposition ($\sqrt{k}$ baby steps, then giant
  steps, then the remaining powers in parallel blocks on the thread
  pool); `nmod_poly_powers_mod_naive` computes them one after the other.

[](nmod_poly/compose_mod) describes these functions, with Rabin's
irreducibility test as an example.

::: functions nmod_poly_compose nmod_poly_compose_series nmod_poly_compose_mod nmod_poly_compose_mod_brent_kung nmod_poly_compose_mod_brent_kung_preinv nmod_poly_precompute_matrix nmod_poly_taylor_shift nmod_poly_revert_series :::
