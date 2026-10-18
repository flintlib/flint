---
title: Rounding in nfloat
summary: Which nfloat operations are exact, which are correctly rounded in every mode or only with directed rounding, which return valid bounds, and which ignore the rounding flags.
order: 2
---

The rounding mode of an `nfloat` context is set with the flags
`NFLOAT_RND_FLOOR` and `NFLOAT_RND_CEIL` of `nfloat_ctx_init`. Without
either flag (the *default mode*), results are rounded in an unspecified
direction for efficiency, usually toward zero. With one of them, the
operations listed below as correctly rounded return the floor or the
ceiling of the exact result, and the others listed as bounded return a
lower or upper bound, which makes it possible to compute rigorous
enclosures with pairs of contexts. Directed rounding is experimental
and may not work correctly together with `NFLOAT_ALLOW_UNDERFLOW`. The
matrix multiplication of `nfloat_mat_mul` falls back to the classical
algorithm in contexts with directed rounding.

## Exact operations {#exact}

These return the exact value in every mode: `nfloat_set`, `nfloat_neg`,
`nfloat_abs`, `nfloat_mul_2exp_si` (barring overflow and underflow),
`nfloat_floor`, `nfloat_ceil`, `nfloat_nint`, `nfloat_trunc`,
`nfloat_sgn`, and the conversions out of `nfloat` that are exact
(`nfloat_get_arf`, `nfloat_get_fmpz`).

## Correctly rounded in every mode {#correct}

In the default mode the result is the exact value truncated toward
zero; with `NFLOAT_RND_FLOOR` or `NFLOAT_RND_CEIL` it is the correctly
rounded floor or ceiling:

- the conversions `nfloat_set_ui`, `nfloat_set_si`, `nfloat_set_fmpz`,
  `nfloat_set_fmpq`, `nfloat_set_d`, `nfloat_set_arf`, and
  `nfloat_set_other` with an `fmpz`, `fmpq`, `arf_t` or `nfloat`
  operand;
- `gr_mul_ui`, `gr_mul_si`, `nfloat_div_ui`, `nfloat_div_si`;
- `nfloat_sqrt`;
- `nfloat_hypot` and `nfloat_complex_abs`.

## Correctly rounded with directed rounding {#directed}

With `NFLOAT_RND_FLOOR` or `NFLOAT_RND_CEIL` these give the correctly
rounded floor or ceiling; in the default mode the error is as noted:

- `nfloat_add`, `nfloat_sub`, and `gr_add_ui`, `gr_add_si`,
  `gr_sub_ui`, `gr_sub_si`: less than 2 ulp (the smaller operand is
  truncated when it is aligned);
- `nfloat_mul`, `nfloat_sqr`: the truncation, except that from three
  words on, where an approximate high product is used, the result can
  rarely be one ulp smaller in magnitude;
- `nfloat_div`, `nfloat_inv`: the truncation, or for divisors of three
  or more nonzero words possibly one ulp larger in magnitude (an error
  of less than 1 ulp);
- `nfloat_rsqrt`: less than $1 + 2^{-37}$ ulp;
- the vector versions of these operations (for example
  `_gr_vec_mul_scalar_ui`), with the exception noted for
  `_nfloat_vec_div_scalar`.

## Valid bounds with directed rounding {#bounds}

With `NFLOAT_RND_FLOOR` or `NFLOAT_RND_CEIL` these return a lower or an
upper bound for the exact value, but not in general the correctly
rounded one:

- `nfloat_addmul`, `nfloat_submul`, `gr_addmul_ui`, `gr_addmul_si`
  (two roundings, so that the result can be far from the exact value,
  measured in ulp of the result, when the terms cancel);
- `_nfloat_vec_dot`, `_nfloat_vec_dot_rev` (the fixed-point sum is
  perturbed by a bound on its error: $n + 2$ units of the last word of
  the $(n + 1)$-word accumulator per term, and 4 per term for $n \le 2$), and matrix multiplication;
- the real elementary functions (`nfloat_exp`, `nfloat_log`,
  `nfloat_sin`, `nfloat_atan`, `nfloat_pow`, the less common functions,
  and so on): within a few ulp at the function precision, and about one
  ulp in the default mode.

## Not respecting the rounding mode {#ignored}

These ignore the rounding flags, and the direction of the error is
arbitrary:

- fused operations and generic algorithms other than those above;
- `nfloat_set_str` (which rounds to nearest);
- `nfloat_set_other` from types other than those listed above;
- mixed operations with an operand of another ring, `fmpz` or `fmpq`
  (`gr_add_other`, `gr_mul_fmpz`, ...);
- the transcendental functions `nfloat_gamma` and `nfloat_zeta`;
- all `nfloat_complex` operations other than `nfloat_complex_abs`; the
  complex elementary functions have errors of a few ulp relative to the
  magnitude of the result in the default mode.

Results that are not correctly rounded are not guaranteed to be
identical on 64-bit and 32-bit machines.
