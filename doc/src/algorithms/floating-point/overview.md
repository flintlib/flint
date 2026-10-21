---
title: Floating-point arithmetic in FLINT
kind: topic
summary: A map of FLINT's floating-point types - arf and arb, nfloat, mp_real, decimal, dfloat and machine doubles - with what each represents, its precision range, whether it is rigorous, and which to use.
order: 0
---

FLINT has several floating-point types, which differ in their radix,
their precision range, their memory layout, and whether they carry
rigorous error bounds. They are not alternatives for the same job so
much as layers: arbitrary-precision balls for rigorous numerics,
fixed-size types for speed in the generic algorithms, and limb- or
double-level kernels on which the others are built.

## The types {#types}

| Type | Module | Radix and precision | Error control | Generic ring |
|---|---|---|---|---|
| `arf_t`, `acf_t` | `arf`, `acf` | binary, any number of bits | none (exact or rounded with a rounding mode, like MPFR) | `gr_ctx_init_real_float_arf` |
| `arb_t`, `acb_t` | `arb`, `acb` | binary, any number of bits | rigorous balls with `mag_t` radius | `gr_ctx_init_real_arb` |
| `nfloat` | `nfloat` | binary, 1 to 66 whole words | none (correctly rounded or within a few ulps; optional directed rounding) | `nfloat_ctx_init` |
| `mp_real_t` | `mp_real` | radix $2^{64}$, whole words | rigorous balls, radius in ulps of the last word | none (a backend) |
| `decfloat_t`, `decball_t` | `decimal` | radix 10, any number of digits | rounded, or rigorous decimal balls | `gr_ctx_init_decfloat`, `gr_ctx_init_decball` |
| `dN_t`, `dNb_t` | `dfloat` | sums of 1 to 4 doubles (53 to 212 bits) | none, or rigorous balls | `gr_ctx_init_dfloat` |
| `double`, `di_t` | `double_extras`, `d_vec`, `d_mat`, `double_interval` | 53 bits | none, or intervals | none |

All of them except the machine doubles have complex counterparts.

## Arbitrary precision: arf and arb {#arb}

`arf_t` is FLINT's general binary floating-point type: the mantissa is
stored inline for up to two words and on the heap beyond, the exponent
is an `fmpz` (unbounded), and every operation takes a precision and a
rounding mode. `arb_t` is a ball with an `arf_t` midpoint and a `mag_t`
radius, and is the basis of FLINT's rigorous numerics (special
functions, polynomial root isolation, integration). These are the types
to use when precision must vary during a computation, exceeds a few
thousand bits, or results must be certified; see
[](algorithms/ball-arithmetic/overview).

## Fixed precision for generic algorithms: nfloat {#nfloat}

`nfloat` numbers have a precision of $n$ whole words, fixed by the
context, and a bounded exponent. Elements are flat arrays of $n + 2$
words with no heap memory, so vectors and matrices are contiguous and
each operation dispatches to a kernel for its size. Over `gr_poly` and
`gr_mat`, approximate root-finding and linear solving run about twice as
fast as over `arf_t` or MPFR [@Joh2025gr]. Results are not certified,
but directed rounding contexts give bounds for the arithmetic, dot
products and elementary functions. See [](algorithms/floating-point/nfloat).

## Limb-level kernels and balls: mp_real {#mp-real}

`mp_real` provides multiple-precision arithmetic in radix $2^{64}$ at
two levels: fixed-point numbers on word arrays, and balls `mp_real_t`
whose radius is counted in units of the last word. It is the backend
for evaluating elementary functions, constants and hypergeometric
series quickly at any precision: the `nfloat` elementary functions run
on its fixed-point kernels. See
[](algorithms/floating-point/mp-real-balls),
[](algorithms/floating-point/elementary-functions),
[](algorithms/floating-point/diophantine-reduction),
[](algorithms/floating-point/hypergeometric-series) and
[](algorithms/floating-point/constants).

## Decimal arithmetic {#decimal}

The `decimal` module provides floating-point numbers and balls in radix
10, with precision measured in decimal digits and mantissas stored in
words of radix $10^{19}$, so that decimal input and output and decimal
rounding are exact. It is the right choice when decimal inputs must be
represented without conversion error or results must be correctly
rounded to a number of decimal digits; binary types are faster
otherwise. See [](algorithms/floating-point/decimal-arithmetic).

## Machine doubles and their extensions: dfloat {#dfloat}

`dfloat` represents numbers as unevaluated sums of 1 to 4 doubles, with
plain formats (approximate, like `double`) and ball formats with a
double radius (rigorous). All operations are built from double
arithmetic and fused multiply-adds, with SIMD for vectors, so that at
53 to 106 bits they are several times faster than `nfloat` and
`arb_t`. The balls track every rounding error exactly, so exact
operations stay exact. They are used for the main sums of the
Riemann–Siegel and Platt methods for the zeta function. See
[](algorithms/floating-point/double-word) and
[](algorithms/floating-point/dfloat-elementary).

The plain `double` helpers (`double_extras`, `d_vec`, `d_mat`) and the
interval type `di_t` of `double_interval` are small support modules,
used for heuristics, tuning and cheap bounds inside other algorithms;
see [](algorithms/floating-point/double-intervals).

## Which to use {#choosing}

| Need | Type |
|---|---|
| Certified results, special functions, any precision | `arb_t`, `acb_t` |
| Certified results at 53 to 212 bits, elementary functions, vectors | `dfloat` balls (`d2b`, `d3b`, ...) |
| Fast approximate linear algebra or root-finding at 64 to 4224 bits | `nfloat` |
| Fast approximate arithmetic at 53 to 212 bits | `dfloat` plain formats |
| Precision varying during the computation, or above 4224 bits | `arf_t` |
| Exact decimal rounding | `decimal` |
| Writing new limb-level numerical kernels | `mp_real` |

Speed comparisons are described in [](dfloat/performance) (which also includes
`nfloat` and `arb_t`); the `nfloat` speedups over `arf_t` are discussed
in [@Joh2025gr].
