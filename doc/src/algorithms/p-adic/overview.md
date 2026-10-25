---
title: p-adic arithmetic in FLINT
kind: topic
summary: The modules for p-adic numbers, their precision models and representations, and which algorithm each function uses for inversion, square roots, the exponential and logarithm, Teichmüller lifts, the Frobenius map, norms and traces.
order: 0
---

FLINT has two independent implementations of $\mathbf{Q}_p$. The
original modules `padic`, `padic_poly`, `padic_mat` and `qadic` work
with any prime $p$, store elements as `fmpz` integers and use a fixed
precision per variable. The newer `padic_radix` works for word-size $p$,
stores the unit in a radix $p^e$ representation, tracks an error term
for every element, and implements the generic ring interface, so that
polynomials, matrices and power series over it are available through
`gr_poly`, `gr_mat` and `gr_series`.

## Modules and representations {#modules}

| Module | Elements | Precision model | $p$ |
|---|---|---|---|
| `padic` | $p^v u$, $u$ an `fmpz` | fixed absolute precision $N$ per variable, inputs exact | any |
| `padic_poly` | $p^v f(X)$, $f \in \Z[X]$ | one $v$ and one $N$ for all coefficients | any |
| `padic_mat` | $p^v M$, $M$ an `fmpz_mat` | one $v$ and one $N$ for all entries | any |
| `qadic` | $p^v g(X) \bmod f$ in $\mathbf{Q}_q = \mathbf{Q}_p[X]/(f)$ | as `padic_poly` | any |
| `padic_radix` | $p^v u + O(p^N)$, $u$ a `radix_integer_t` | ball with its own $N$; context caps on absolute and relative precision | word-size |

In the fixed-precision modules every function computes the exact result
for the stored inputs and reduces it modulo $p^N$ for the output; no
loss of precision is detected. In `padic_radix` every function applies
the precision rules of [](algorithms/p-adic/precision#propagation) to
the error terms, and exact elements remain exact when the result fits
the context precision. Both models are described in
[](algorithms/p-adic/precision). Unramified extensions exist only in
the fixed-precision model (`qadic`); ramified extensions are not
implemented.

The storage differs as well. An `fmpz` unit needs
$\lceil (N - v)\log_2 p \rceil$ bits, and reduction modulo $p^{N - v}$
is a division by an integer. A `radix_integer_t` stores the base-$p^e$
digits, with $p^e$ chosen as large as possible (for example $7^{22}$
on a 64-bit machine), so
reduction modulo a power of $p^e$ is a truncation, valuations are found
by inspecting the lowest digit, and no conversion from binary is
needed for printing.

## Algorithms {#algorithms}

| Operation | Method | Functions | Details |
|---|---|---|---|
| Inverse of a unit | extended Euclid modulo $p$, Newton lifting $z \gets 2z - uz^2$ | `padic_inv`, `_padic_inv`, `qadic_inv`, `padic_radix_inv` | [](algorithms/p-adic/hensel-lifting#inverse) |
| Square root | root modulo $p$ (or Artin–Schreier modulo 8 for $p = 2$ in $\mathbf{Z}_q$), lifting of the inverse square root | `padic_sqrt`, `qadic_sqrt`, `padic_radix_sqrt`, `padic_radix_rsqrt` | [](algorithms/p-adic/hensel-lifting#sqrt) |
| Exponential | rectangular splitting, or bit-burst splitting with binary splitting | `padic_exp*`, `qadic_exp*`, `padic_radix_exp*` | [](algorithms/p-adic/exp-log) |
| Logarithm | rectangular splitting, product decomposition with binary splitting, or the Satoh–Skjernaa–Taguchi reduction | `padic_log*`, `qadic_log*`, `padic_radix_log*` | [](algorithms/p-adic/exp-log) |
| Teichmüller lift | Newton iteration for $t^q - t$ | `padic_teichmuller`, `qadic_teichmuller` | [](algorithms/p-adic/hensel-lifting#teichmuller) |
| Frobenius $\Sigma^e$ | Newton lifting of $\Sigma^e(X)$, composition by rectangular splitting | `qadic_frobenius` | [](algorithms/p-adic/frobenius) |
| Trace in $\mathbf{Q}_q$ | power sums by Newton's identities | `qadic_trace` | [](algorithms/p-adic/unramified-extensions#trace) |
| Norm in $\mathbf{Q}_q$ | Sylvester determinant, or $\exp(\operatorname{Tr}\log x)$ | `qadic_norm`, `qadic_norm_resultant`, `qadic_norm_analytic` | [](algorithms/p-adic/unramified-extensions#norm) |
| Products of polynomials and matrices | exact products over $\Z$, then reduction | `padic_poly_mul`, `padic_mat_mul`, `qadic_mul` | [](algorithms/polynomials/multiplication-integers) |
| Dot products, matrix products | delayed carries over radix limbs; limb-aligned packed kernel, or products modulo $B^S$ (`nmod_mat`, `mpn_mod_mat`, `fmpz_mat`); precisions by min-plus products | `padic_radix_dot`, `padic_radix_mat_mul` | [](algorithms/p-adic/block-products) |
| Polynomial products over `padic_radix` | packed kernel, Kronecker substitution in radix $p^e$, or products modulo $B^S$ (`nmod_poly`, `mpn_mod_poly`); precisions by exact min-plus convolutions | `_padic_radix_poly_mulmid` | [](algorithms/p-adic/block-products#polynomials) |
| Determinant over `padic_radix` | full pivoting on valuations, blocked by valuation levels (optimal precision) | `padic_radix_mat_det` | [](algorithms/p-adic/linear-algebra#det) |
| Linear systems, inverse over `padic_radix` | Gauss–Jordan elimination by valuation levels (optimal precision) | `padic_radix_mat_nonsingular_solve`, `padic_radix_mat_inv` | [](algorithms/p-adic/linear-algebra#solve) |
| Characteristic polynomial over `padic_radix` | Berkowitz for cyclic matrices, otherwise an exact lift with the precision of the differential | `padic_radix_mat_charpoly` | [](algorithms/p-adic/linear-algebra#charpoly) |

`padic_poly_mul` and `padic_mat_mul` multiply the integer parts with
`fmpz_poly_mul` and `fmpz_mat_mul` and reduce once; the valuations
simply add. `padic_poly_inv_series` lifts a power series inverse by
Newton iteration over $\Z/p^{N}\Z$. Polynomial and matrix algorithms not
listed here (division, GCD, factorisation, linear solving) are not
implemented in the fixed-precision modules; over `padic_radix` the
generic algorithms of `gr_poly` and `gr_mat` apply, with the precision
loss tracked by the ball arithmetic. For products, determinants, linear
systems and characteristic polynomials `padic_radix` has its own
methods, which give the optimal precision
([](algorithms/p-adic/linear-algebra)).

## Choosing a module {#choosing}

- For new code with word-size $p$, `padic_radix` is faster and reports
  the precision of every result; it is the only choice if precision
  loss must be detected, and the only module with polynomials and
  matrices that support division and linear algebra (through the generic
  rings).
- For multiword primes, for unramified extensions, and for code that
  needs `padic_poly`, `padic_mat` or the Frobenius map, the
  fixed-precision modules are needed. Their results are exact functions
  of the stored inputs, so the caller has to bound the precision loss
  by the rules of [](algorithms/p-adic/precision).
