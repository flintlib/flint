---
title: Zech logarithms
summary: In a small field $\mathbb{F}_q$ with a primitive element $g$, storing each element as its discrete logarithm turns multiplication into an addition modulo $q - 1$, and a table of the Zech logarithms $Z(k)$ with $1 + g^k = g^{Z(k)}$ turns addition into one lookup; the table is built in $O(q)$ operations, which also tests whether $g$ is primitive.
complexity: $O(1)$ word operations for every field operation; $O(q)$ field operations and $2q + p$ words of memory to build the tables.
order: 20
---

The multiplicative group $\mathbb{F}_q^\times$ of a finite field is
cyclic of order $q - 1$. If $g$ is a generator, every nonzero element is
$a = g^k$ for a unique $k \in [0, q - 2]$, its discrete logarithm, and
multiplication, division and powering of elements become addition,
subtraction and multiplication of logarithms modulo $q - 1$. Addition
needs a table: Jacobi introduced, and Zech tabulated, the function $Z$
with

$$
1 + g^k = g^{Z(k)},
$$

now called the Zech (or Jacobi) logarithm [@LN1997] [@Hub1990].
The representation is attractive whenever a table of size $q$ is
affordable; `fq_default` uses it for $q < 2^{16}$, and the module
`fq_zech` allows $q < 2^{\bits}$.

## Arithmetic with logarithms {#arithmetic}

In `fq_zech` an element is stored as one word: its logarithm $k$, or
$q - 1$ for zero (playing the role of $\log 0 = \infty$). Then, with all
logarithms taken modulo $q - 1$ and the special value for zero handled
separately,

- $g^k g^l = g^{k + l}$, $(g^k)^{-1} = g^{-k}$, $(g^k)^e = g^{ke}$;
- the Frobenius map is $\sigma^e(g^k) = g^{k p^e}$ and the $p$-th root
  is $g^{k p^{d-1}}$;
- $-1 = g^{(q-1)/2}$ for odd $p$ (the unique element of order 2), and
  $-1 = 1 = g^0$ for $p = 2$, so negation adds $(q - 1)/2$ or nothing;
- for odd $p$, $g^k$ is a square if and only if $k$ is even, since the
  squares form the subgroup of index 2; its square roots are
  $\pm g^{k/2}$. For $p = 2$ every element is a square, and the square
  root of $g^k$ is $g^{k/2}$ or $g^{(k + q - 1)/2}$ ($q - 1$ is odd);
- for $i \ge j$, $g^i + g^j = g^j (1 + g^{i - j}) = g^{j + Z(i - j)}$,
  where $Z(i - j)$ is "infinite" (the sum is zero) when $g^{i-j} = -1$;
  subtraction is addition of $g^{j + (q-1)/2}$.

So every operation except table construction costs a few word operations
and one memory access. When the table does not fit in the cache, the
memory access dominates: for $q \approx 10^6$ the polynomial
representation `fq_nmod` is already faster (see
[](algorithms/finite-fields/overview)).

## Building the tables {#tables}

Let $f$ be the defining polynomial and $g$ the class of $x$ in
$\mathbb{F}_p[x]/(f)$. To index elements by integers, write an element
in the polynomial representation as $\sum_{j < d} c_j g^j$ with
$c_j \in [0, p)$ and map it to the integer
$\nu(a) = \sum_j c_j p^j \in [0, q)$, its coefficients read as digits in
base $p$. This is a bijection between $\mathbb{F}_q$ and $[0, q)$, with
$\nu(0) = 0$ and $\nu(1) = 1$.

::: algorithm Zech tables {#alg-tables}
**Input:** a monic irreducible $f$ of degree $d$ over $\mathbb{F}_p$.

**Output:** the tables $E[k] = \nu(g^k)$ and $Z(k)$ for $0 \le k < q - 1$,
and $L[c] = \log c$ for $c \in \mathbb{F}_p^\times$; or failure if $g$ is
not primitive.

1. Set $T[n] \gets \infty$ for $0 \le n < q$, and $r \gets 1$.
2. For $k = 0, 1, \ldots, q - 2$: let $n = \nu(r)$; if $T[n] \ne \infty$,
   fail. Set $T[n] \gets k$ and $E[k] \gets n$; if $r$ is a constant
   $c$, set $L[c] \gets k$. Set $r \gets r g$.
3. For each $n \in [0, q)$: let $n' = n + 1 - p$ if $n \equiv p - 1
   \pmod p$ and $n' = n + 1$ otherwise; set $Z(T[n]) \gets T[n']$ (with
   $T[0] = \infty$ for the element 0).
:::

::: theorem {#thm-tables}
[](#alg-tables) fails if and only if $g$ is not a primitive element of
$\mathbb{F}_q$. Otherwise it computes the tables correctly, with $q - 1$
multiplications in $\mathbb{F}_q$ and $O(q)$ further word operations.
:::

::: proof
The loop in step 2 visits $1, g, g^2, \ldots, g^{q-2}$. These are pairwise
distinct if and only if the order of $g$ is at least $q - 1$, that is, if
$g$ is primitive; since $\nu$ is injective and $g^k \ne 0$, a repeated
value of $\nu$ occurs exactly when two of the powers coincide. If $g$ is
primitive, every nonzero element is reached once, so $T[\nu(a)] = \log a$
for every $a \ne 0$, and $E$ and $L$ are the stated tables. In step 3,
$n'$ is $\nu(a + 1)$ for the element $a$ with $\nu(a) = n$: adding 1
increases the constant coefficient $c_0 = n \bmod p$ by one, except that
$p - 1$ wraps to 0, which subtracts $p - 1$ from $n$. Hence for
$a = g^k$, $T[n'] = \log(1 + g^k) = Z(k)$, with the value $\infty$ when
$1 + g^k = 0$.
:::

The tables take $q$ words for $Z$, $q$ for $E$ and $p$ for $L$; the
temporary table $T$ takes another $q$ words. Because the construction
fails for a non-primitive $f$, it doubles as a primitivity test:
`fq_zech_ctx_init_fq_nmod_ctx_check` returns 0 in that case. Conway
polynomials are primitive by definition, and other `fq_zech` contexts are
created with random primitive polynomials
([](algorithms/finite-fields/defining-polynomials)).

## Conversions, norm and trace {#conversions}

The table $E$ gives the polynomial representation of $g^k$ directly: its
coefficients are the base-$p$ digits of $E[k]$, which is how
`fq_zech_get_fq_nmod` and the polynomial multiplication of `fq_zech_poly`
convert elements. Conversely, $\sum_j c_j g^j$ is computed as a sum of
the terms $g^{L[c_j] + j}$ with $d - 1$ Zech additions.

The norm of $g$ is $N(g) = \prod_{i < d} g^{p^i} = g^{(q-1)/(p-1)}$, and
it is $(-1)^d f(0)$ for the monic minimal polynomial $f$ of $g$; it is a
generator of $\mathbb{F}_p^\times$, since $g^{(q-1)/(p-1)}$ has order
$p - 1$. The context stores it (`prime_root`), and then:

- $N(g^k) = N(g)^k$, one exponentiation modulo $p$ (`fq_zech_norm`);
- the elements of $\mathbb{F}_p^\times$ are the powers $g^{m(q-1)/(p-1)}
  = N(g)^m$, so an element of the prime field with logarithm $l$
  corresponds to the integer $N(g)^{l (p-1)/(q-1)} \bmod p$;
  `fq_zech_trace` sums the $d$ conjugates $g^{k p^i}$ with Zech additions
  and converts the result in this way.

## Use in FLINT {#use}

- The tables are built by `fq_zech_ctx_init_fq_nmod_ctx_check`, which is
  called by all `fq_zech_ctx_init` functions and by `fq_default` when it
  chooses the Zech representation.
- `fq_zech_add`, `fq_zech_sub`, `fq_zech_neg`, `fq_zech_mul`,
  `fq_zech_inv`, `fq_zech_pow`, `fq_zech_frobenius`, `fq_zech_pth_root`,
  `fq_zech_sqrt`, `fq_zech_is_square`, `fq_zech_norm` and `fq_zech_trace`
  implement the formulas above.
- The generic `gr` context of `fq_zech` has a dot product working
  directly with logarithms, used by the classical polynomial
  multiplication of `fq_zech_poly` (through `_gr_poly_mullow_classical`);
  long polynomial products convert to the polynomial representation
  ([](fq_poly/mul)).
