---
title: fmpz_extras.h
header: fmpz_extras.h
description: Inline fast paths for fmpz arithmetic on small values, and helpers for exponents, limb access and memory accounting.
category: Integers
order: 3
---

This module provides utilities for [`fmpz`](fmpz) integers that are
mainly used inside FLINT, in particular by the ball arithmetic modules
(`mag`, `arf`, `arb`), whose exponents are `fmpz_t` values that are
nearly always small. Most functions are `static inline` and put the fast
path for inline operands directly into the caller's code, falling back
to the corresponding `fmpz` function otherwise; the guide
[](fmpz/representation#underscore) explains the representation they rely
on. The remaining functions round, compare and convert in ways the
`fmpz` module does not offer.

## Memory {#memory}

::: function
slong fmpz_allocated_bytes(const fmpz_t x)

Returns the number of bytes allocated on the heap for $x$: 0 for an inline
value, and the size of the `mpz` structure plus its allocated limbs
otherwise. Add `sizeof(fmpz)` for the total size of the object.
:::

## Inline arithmetic {#inline}

These functions give the same results as the `fmpz` functions they
replace, for all inputs. If the operands are inline (and the word-size
operand lies in $[$`COEFF_MIN`, `COEFF_MAX`$]$), the result is computed
with one machine addition; otherwise the out-of-line `fmpz` function is
called. Inlining enlarges the code and
slightly slows down the general case, so these functions should be used
only where operands are almost always small.

::: function
void fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)
void fmpz_add_si_inline(fmpz_t z, const fmpz_t x, slong y)
void fmpz_add_ui_inline(fmpz_t z, const fmpz_t x, ulong y)

Set $z$ to $x + y$.
:::

::: function
void fmpz_sub_si_inline(fmpz_t z, const fmpz_t x, slong y)

Sets $z$ to $x - y$.
:::

::: function
void fmpz_add2_fmpz_si_inline(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)

Sets $z$ to $x + y + c$.
:::

::: function
slong _fmpz_sub_small(const fmpz_t x, const fmpz_t y)

Returns $x - y$ as an `slong` when it can be represented cheaply: exactly
if $x$ and $y$ are both inline or if $|x - y| \le$ `COEFF_MAX`, and
otherwise `WORD_MAX` or `-WORD_MAX` with the sign of $x - y$. The sign of
the result is always correct, which suffices to compare exponents and to
decide whether two floating-point numbers overlap.
:::

::: function
slong _fmpz_sub_small_large(const fmpz_t x, const fmpz_t y)

The out-of-line part of `_fmpz_sub_small`, for operands that are not both
inline.
:::

::: function
void _fmpz_set_si_small(fmpz_t x, slong v)

Sets $x$ to $v$, which must satisfy `COEFF_MIN` $\le v \le$ `COEFF_MAX`
(not checked), so that no promotion is needed.
:::

::: function
slong _fmpz_size(const fmpz_t x)

Returns the number of limbs of $|x|$, counting an inline value
(including zero) as one limb. Unlike `fmpz_size`, it returns 1 for 0.
:::

## Convenience functions {#convenience}

::: function
void fmpz_adiv_q_2exp(fmpz_t z, const fmpz_t x, flint_bitcnt_t exp)

Sets $z$ to $x / 2^{\mathit{exp}}$ rounded away from zero.
:::

::: function
void fmpz_ui_mul_ui(fmpz_t x, ulong a, ulong b)

Sets $x$ to the product $ab$ of two words, which may need two limbs.
:::

::: function
void fmpz_max(fmpz_t z, const fmpz_t x, const fmpz_t y)
void fmpz_min(fmpz_t z, const fmpz_t x, const fmpz_t y)

Set $z$ to the maximum, respectively the minimum, of $x$ and $y$.
:::

## Low-level conversions {#conversions}

::: macro
FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, zv)

Given the `fmpz` value `zv` (a word, not an `fmpz_t`), sets `zptr` to a
pointer to the limbs of $|zv|$, `zn` to their number and `zsign` to 1 if
`zv` is negative and 0 otherwise. For an inline value the absolute value
is stored in the `ulong` variable `ztmp`, which `zptr` then points to, so
`ztmp` must stay in scope while `zptr` is used. Zero is not handled
specially: it gives `zn` $= 1$ and a zero limb.
:::

::: function
void fmpz_lshift_mpn(fmpz_t z, nn_srcptr src, slong n, int negative, flint_bitcnt_t shift)

Sets $z$ to the integer given by the $n \ge 1$ limbs `src` (with a
nonzero top limb), negated if `negative` is 1, and multiplied by
$2^{\mathit{shift}}$. Used to convert mantissas of floating-point
numbers to integers. The result is demoted if it is small.
:::
