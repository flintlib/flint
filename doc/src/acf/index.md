---
title: acf.h
header: acf.h
description: Complex floating-point numbers as pairs of arf_t, for approximate computation
category: Real and complex numbers
order: 5
---

An `acf_t` is a complex floating-point number $a + bi$ whose real and
imaginary parts are `arf_t` numbers, without error bounds. It is the
complex counterpart of `arf_t` and is used where complex arithmetic is
needed but error bounds are not: in approximate (non-rigorous)
algorithms such as the floating-point stage of eigenvalue or root
computations that are certified afterwards, and as the element type of the
generic ring of complex floating-point numbers
(`gr_ctx_init_complex_float_acf`).

Addition, subtraction and multiplication are
correctly rounded in each part, like the `arf` operations; division,
square root and dot products are only approximate. All functions allow
aliasing of inputs and outputs.

## Types {#types}

::: type
acf_struct
acf_t
acf_ptr
acf_srcptr

An `acf_struct` has two `arf_struct` fields, the real and imaginary
parts. `acf_t` is an array of one `acf_struct`, passed by reference;
`acf_ptr` and `acf_srcptr` point to arrays.
:::

::: macro
acf_realref(x)
acf_imagref(x)

Pointers to the real and imaginary parts of `x` as `arf_t`.
:::

::: function
arf_ptr acf_real_ptr(acf_t z)
arf_ptr acf_imag_ptr(acf_t z)

Function versions of `acf_realref` and `acf_imagref`, for language
bindings.
:::

## Memory management {#memory}

::: function
void acf_init(acf_t x)
void acf_clear(acf_t x)

Initialise `x` to zero, and free the memory used by `x`.
:::

::: function
acf_ptr _acf_vec_init(slong n)
void _acf_vec_clear(acf_ptr v, slong n)

Allocate a vector of `n` numbers initialised to zero, and free such a
vector.
:::

::: function
void acf_swap(acf_t z, acf_t x)

Swaps `z` and `x` without copying.
:::

::: function
slong acf_allocated_bytes(const acf_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(acf_struct)`.
:::

## Basic operations {#basic}

::: function
void acf_set(acf_t z, const acf_t x)
void acf_neg(acf_t z, const acf_t x)

Set `z` to `x`, respectively $-x$, exactly.
:::

::: function
int acf_set_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_neg_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to `x`, respectively $-x$, with both parts rounded to `prec`
bits in the direction `rnd`. Bit 0 of the return value is the inexact
flag of the real part and bit 1 that of the imaginary part.
:::

::: function
int acf_equal(const acf_t x, const acf_t y)

Returns whether `x` and `y` are equal (both parts equal in the sense of
`arf_equal`).
:::

::: function
slong acf_bits(const acf_t x)

Returns the larger of `arf_bits` of the two parts.
:::

::: function
void acf_get_mag(mag_t res, const acf_t x)

Sets `res` to an upper bound for $|x|$.
:::

::: function
void acf_printd(const acf_t x, slong n)

Prints `x` as `a + b*I` with `n` significant decimal digits for each part
(see `arf_printd`).
:::

::: function
void acf_randtest(acf_t x, flint_rand_t state, slong bits, slong mag_bits)

Sets both parts of `x` to random numbers as `arf_randtest` does.
:::

## Arithmetic {#arithmetic}

::: function
int acf_add(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_sub(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_mul(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to $x + y$, $x - y$, respectively $xy$, with each part correctly
rounded to `prec` bits in the direction `rnd`. Bit 0 of the return value
is the inexact flag of the real part and bit 1 that of the imaginary part.

::: details
`acf_mul` uses `arf_complex_mul`, or `arf_complex_sqr` when `x` and `y`
are the same object; each part of the result is the exact $ac - bd$ or
$ad + bc$ rounded once.
:::
:::

## Approximate arithmetic {#approx}

The following functions are not correctly rounded: `rnd` is the direction
of the final rounding, but the intermediate roundings are not specified
and no error bound is given.

::: function
void acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
void acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to an approximation of $1/x$, respectively $x/y$. A real or
purely imaginary divisor is handled by real divisions of the parts;
otherwise $1/(a + bi) = (a - bi)/(a^2 + b^2)$ is used, with $a^2 + b^2$
computed by `arf_sosq`.
:::

::: function
void acf_approx_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of the principal square root of `x`.
:::

::: function
void acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of $s + (-1)^{\mathit{subtract}}
\sum_{i=0}^{\mathit{len}-1} x_i y_i$, where $s$ is `initial` (0 if
`NULL`), with the conventions of `arb_dot`, using a fixed-point
accumulation like that of `acb_approx_dot`.
:::
