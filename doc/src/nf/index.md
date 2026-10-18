---
title: nf.h
header: nf.h
description: Number fields $\mathbf{Q}[x]/(f)$ given by a defining polynomial over $\mathbf{Q}$, with precomputed data for the arithmetic of `nf_elem.h`.
category: Algebraic numbers and number fields
order: 2
---

This module defines the type `nf_t` of a number field
$K = \mathbf{Q}(\theta) = \mathbf{Q}[x]/(f)$, given by a defining
polynomial $f \in \mathbf{Q}[x]$ of degree $d \ge 1$, and computes the
data that the element arithmetic of `nf_elem.h` needs: the
representation used for elements (special ones for degree 1 and 2), a
precomputed inverse of the leading coefficient, the powers
$x^i \bmod f$ used for reduction, and the power sums of the roots of $f$
used for traces. The module goes back to the ANTIC library
[@Har2015]. The generic ring interface for number fields is
`gr_ctx_init_nf` in `gr.h`; algebraic number fields as subfields of
$\mathbf{C}$ are handled by `qqbar.h` and `ca.h`.

::: summary algorithms/number-theory/number-fields :::

## Types {#types}

::: type
nf_struct
nf_t

A number field. The fields are `pol` (the defining polynomial, an
`fmpq_poly_t`), `pinv` (a precomputed inverse of the leading coefficient
of the numerator of `pol`, used for division in the non-monic generic
case), `powers` (precomputed reductions of $x^{d}, \ldots, x^{2d-2}$
modulo `pol` for $3 \le d < $ `NF_POWERS_CUTOFF`), `traces` (the power
sums $S_k = \sum_i \theta_i^k$ of the roots $\theta_i$ of `pol` for
$0 \le k < d$, as an `fmpq_poly_t`) and `flag`.
:::

::: macro
NF_GENERIC
NF_MONIC
NF_LINEAR
NF_QUADRATIC
NF_GAUSSIAN

Bits of the field `flag`. `NF_MONIC` is set if `pol` has integer
coefficients and leading coefficient 1, `NF_LINEAR` for $d = 1$,
`NF_QUADRATIC` for $d = 2$, and `NF_GAUSSIAN` in addition for
$f = x^2 + 1$ exactly; `NF_GENERIC` is the value 0. Elements of linear
and quadratic fields use the representations `lnf_elem_t` and
`qnf_elem_t` of `nf_elem.h`.
:::

::: macro
NF_POWERS_CUTOFF

The bound 30 on the length of the defining polynomial up to which the
powers of $x$ modulo $f$ are precomputed for reduction.
:::

## Functions {#functions}

::: function
void nf_init(nf_t nf, const fmpq_poly_t pol)

Initialises `nf` for the defining polynomial `pol`, which must have
degree at least 1 (otherwise the function throws). The polynomial is
not checked to be irreducible: with a reducible polynomial the
arithmetic is that of the ring $\mathbf{Q}[x]/(f)$, and the result of
inverting a zero divisor is meaningless (`_nf_elem_invertible_check`
tests for this).

::: details
The power sums are computed by the recurrence of Newton's identities
for a non-monic polynomial [@Coh1996, p. 163], with a common denominator
that is a power of the leading coefficient.
:::
:::

::: function
void nf_init_randtest(nf_t nf, flint_rand_t state, slong len, flint_bitcnt_t bits_in)

Initialises `nf` with a random defining polynomial of length at most
`len` (that is, degree at most $len - 1$) with coefficients of at most
`bits_in` bits, for testing. Linear and quadratic fields, monic
polynomials and $x^2 + 1$ are chosen with increased probability. The
polynomial need not be irreducible. Throws if $len < 2$ or
`bits_in` $< 1$.
:::

::: function
void nf_clear(nf_t nf)

Clears `nf`.
:::

::: function
void nf_print(const nf_t nf)

Prints a description of `nf`, with the defining polynomial in the
variable `x`, to `stdout`, without a final newline.
:::
