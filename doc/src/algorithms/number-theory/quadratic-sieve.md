---
title: The quadratic sieve
summary: Factoring $N$ by finding a congruence of squares from many smooth values of quadratic polynomials, with self-initialising polynomials, a block and bucket sieve, the large prime variation and block Lanczos linear algebra.
complexity: Heuristically $L_N[1/2, 1] = \exp\bigl((1 + o(1)) \sqrt{\log N \log \log N}\bigr)$.
order: 32
---

The quadratic sieve of Pomerance [@Pom1985] is the fastest known method
for factoring integers without small factors of up to about 100 digits.
FLINT implements the *self-initialising multiple polynomial* variant
(SIQS) [@AlfPom1995; @Con1997] with the large prime variation, in the
`qsieve` module; `fmpz_factor` uses it for the cofactors left by trial
division and the [elliptic curve method](algorithms/number-theory/ecm).
The implementation is due to William Hart and Nitin Kumar (2015–2016),
with the block Lanczos code of Jason Papadopoulos (from msieve).

## Congruences of squares {#squares}

::: lemma {#lem-squares}
Let $N$ be odd with at least two distinct prime factors, and let
$X^2 \equiv Y^2 \pmod N$. If $X \not\equiv \pm Y \pmod N$, then
$\gcd(X - Y, N)$ is a nontrivial factor of $N$. If $Y$ is a random square
root of $X^2$ modulo $N$ (with $\gcd(X, N) = 1$), this happens with
probability at least $1/2$.
:::

::: proof
$N$ divides $(X - Y)(X + Y)$ but neither factor, so the GCD is proper.
If $N$ has $r \ge 2$ distinct prime factors, $X^2$ has $2^r$ square roots
modulo $N$ (by the Chinese remainder theorem, as $N$ is odd), of which
only two are $\pm X$.
:::

Dixon's and Kraitchik's idea is to produce such congruences from
*relations*: values $x$ for which $Q(x) = x^2 - N$ factors completely
over a *factor base* $\mathcal B$ of small primes $p_1, \ldots, p_F$ (and
$-1$). If $Q(x_i) = \pm \prod_j p_j^{e_{ij}}$, a subset $S$ of the
relations with $\sum_{i \in S} e_{ij}$ even for all $j$ (and an even
number of negative values) gives $X = \prod_{i \in S} x_i$ and
$Y = \prod_j p_j^{\frac12 \sum_{i \in S} e_{ij}}$ with $X^2 \equiv Y^2$.
Such subsets are the vectors of the null space of the exponent matrix
modulo 2, and $F + 1$ relations guarantee a nonzero one; with $F + 64$
relations one expects 64 independent dependencies, each splitting $N$
with probability at least $1/2$.

## Sieving {#sieving}

Only primes $p$ for which $N$ is a quadratic residue can divide some
$Q(x)$, so $\mathcal B$ consists of those primes (about half of all
primes up to the bound). For such $p$, $p \mid Q(x)$ if and only if
$x \equiv \pm s_p \pmod p$ where $s_p^2 \equiv N$. So instead of trial
dividing every $Q(x)$, one initialises an array indexed by $x$ in an
interval and, for every $p \in \mathcal B$, adds $\log p$ at the
positions $x \equiv \pm s_p \pmod p$. Positions where the accumulated
logarithm is close to $\log |Q(x)|$ are likely to be smooth and are
examined by trial division. Sieving an interval of length $M$ costs
about $M \sum_{p \le b} 2/p \approx 2M \log \log b$ additions, an
amortised constant per position.

A *multiplier* $k$ (squarefree, at most 47) is chosen so that $kN$ is a
square modulo many small primes and $kN \equiv 1 \pmod 8$ if possible,
and $kN$ is used instead of $N$. FLINT chooses it by the Knuth–Schroeppel
function [@Sil1987]: $f(k) = -\frac12 \log k + \sum_p g(p, kN) \log p$,
where $g$ is the expected contribution of $p$ to a sieve value
(`qsieve_knuth_schroeppel`).

## Self-initialising polynomials {#siqs}

The values of $x^2 - kN$ grow linearly with $x$, so the yield drops
along the interval. The multiple polynomial variant [@Sil1987] uses many
polynomials $Q(x) = (Ax + B)^2 - kN$ with $B^2 \equiv kN \pmod A$, for
which $Q(x) = A (Ax^2 + 2Bx + C)$ with $C = (B^2 - kN)/A$. Sieving
$g(x) = Ax^2 + 2Bx + C$ over $[-M, M]$, the values are at most about
$M \sqrt{kN/2}$ if $A \approx \sqrt{2kN}/M$, much smaller than with a
single polynomial.

The self-initialising variant [@AlfPom1995] chooses
$A = q_1 \cdots q_s$ as a product of $s$ factor base primes. Then
$B^2 \equiv kN \pmod A$ has $2^s$ solutions, $2^{s-1}$ up to sign, each
giving a polynomial. They are $B = \sum_{\nu} \pm B_\nu$ with
$B_\nu = (A/q_\nu) \gamma_\nu$, $\gamma_\nu \equiv \sqrt{kN}
(A/q_\nu)^{-1} \pmod {q_\nu}$, and enumerating the signs in Gray code
order changes $B$ by $\pm 2 B_\nu$ at each step.

::: lemma {#lem-roots}
Let $p \in \mathcal B$ with $p \nmid A$. The roots of $g$ modulo $p$ are
$x \equiv A^{-1}(\pm s_p - B) \pmod p$. If $B$ changes to $B \mp 2B_\nu$,
the roots change by $\pm 2 B_\nu A^{-1} \bmod p$.
:::

::: proof
$A g(x) = (Ax + B)^2 - kN \equiv 0 \pmod p$ iff $Ax + B \equiv \pm s_p$,
and $A$ is invertible modulo $p$. The second statement follows by
substituting the new $B$.
:::

So after precomputing $A^{-1} \bmod p$ and $2 B_\nu A^{-1} \bmod p$ for all
$p$ and $\nu$ (the *initialisation* for $A$), each new polynomial costs
one addition per factor base prime. FLINT picks the $q_\nu$ in a range
of the factor base chosen so that $s$ primes of that size give $A$ within
a factor of 2 of its target (following msieve); for $s > 3$ the first
$s - 1$ are chosen at random and the last by binary search so that the
product is close to the target (`qsieve_next_A`).

## The sieve in FLINT {#sieve}

The interval $[-M, M)$ has `sieve_size` bytes, each holding a scaled
logarithm, processed in blocks of $2^{19}$ bytes that fit in the L2
cache. The primes are treated according to their size:

- the smallest primes (the first `small_primes`, typically 4 to 30) are
  not sieved at all; their contribution is compensated by the threshold;
- primes below 4096, which hit a block many times, are sieved over
  sub-blocks of $2^{15}$ bytes that fit in the L1 cache;
- primes up to the block size are sieved directly over the block;
- larger primes hit each block at most a few times. They are processed
  once per polynomial by a *bucket sieve*: every hit is appended to a
  bucket for the sub-block it falls in, and the buckets are applied to the
  sub-blocks in turn. The buckets are kept, so that when a candidate is
  examined, its large factor base primes are read off the buckets instead
  of being found by trial division.

The sieve values are biased so that a position is a candidate exactly
when the top bit of its byte is set, and the scan tests eight positions
at a time. Candidates are factored by trial division with
word-size multiplications by precomputed inverses (in batches of 64
primes). The polynomials for one $A$ are distributed over the threads.

## Large primes and cycles {#large-primes}

A candidate that factors over $\mathcal B$ except for one prime $\ell$
below a bound $L$ (60 times the largest factor base prime) is kept as a
*partial relation*. Two partial relations with the same $\ell$ multiply
to a relation in which $\ell$ occurs squared, which is as good as a full
relation. More generally, consider the graph whose vertices are 1 and the
large primes, with an edge $\{1, \ell\}$ for each partial relation (an
edge $\{\ell_1, \ell_2\}$ for a relation with two large primes).

::: lemma {#lem-cycles}
Every cycle of this graph gives a full relation, and the number of
independent cycles of a graph with $V$ vertices, $E$ edges and $C$
connected components is $E - V + C$.
:::

::: proof
Multiplying the relations along a cycle, every vertex on it is the
endpoint of two of its edges, so every large prime occurs to an even
power and can be moved to $Y$ (and the vertex 1 stands for no prime). The
cycle space of a graph has dimension $E - V + C$, the number of edges
not in a spanning forest.
:::

FLINT counts the cycles as relations arrive with a union-find structure
(`qsieve_add_partial_edge`), so that it knows when there are enough
relations, and when processing them builds a spanning forest and turns
every non-tree edge into a full relation by following the tree paths
(`qsieve_process_relation`). The relation store
and the graph support two large primes per relation, but the sieve
currently produces only single large prime relations.

## Linear algebra and square root {#linalg}

When the full relations and cycles exceed the number of factor base
primes by 64, the relations are put into a sparse matrix over
$\mathbb F_2$, lightly filtered (columns with a singleton row are deleted
repeatedly, since they cannot be part of a dependency), and Montgomery's
block Lanczos algorithm [@Mon1995] finds up to 64 vectors of its null
space in time $O(F \cdot w)$ for a matrix with $w$ nonzero entries, as a
sequence of sparse matrix products with blocks of 64 vectors. For each
dependency, the square root step forms $X = \prod Y_i$ and $Y$ from the
halved exponents modulo $kN$ and takes $\gcd(X - Y, N)$; this stops as
soon as the factors found split $N$ completely into probable primes. If
no dependency splits $N$, the factor base is enlarged and sieving
continues with the relations already found.

## Complexity and parameters {#complexity}

With a factor base bound $b = L_N[1/2, 1/2]$, the heuristic probability
that a value of size $N^{1/2 + o(1)}$ is $b$-smooth is $L_N[1/2, -1/2]$,
so $L_N[1/2, 1/2]$ relations are found by sieving $L_N[1/2, 1]$
positions. The linear algebra takes $L_N[1/2, 1]$ with a sparse method
such as block Lanczos (and $L_N[1/2, 3/2]$ with Gaussian elimination),
so the total is $L_N[1/2, 1]$ [@CraPom2005, §6.1]; in practice the sieve
dominates.

The parameters (number of primes considered for the multiplier, size of
the factor base, number of unsieved small primes, size of the interval and
threshold) come from a table indexed by the bit size of $N$, tuned
experimentally (`qsieve_tune`). For example, for 160-bit numbers (49
digits) the factor base has 1700 primes, for 200 bits (61 digits) 4500
and for 300 bits (91 digits) 140000.

## Use in FLINT {#use}

`qsieve_factor` (and `qsieve_factor_with_tune`) implement the method;
`fmpz_factor` calls it for composite cofactors without small factors,
after checking that they are not perfect powers (which the sieve cannot
split). For timings, see [](fmpz_factor/factor).

::: functions qsieve_factor qsieve_knuth_schroeppel qsieve_process_relation block_lanczos :::
