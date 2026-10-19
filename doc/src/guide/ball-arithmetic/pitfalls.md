---
title: Technical conventions
summary: The conventions specific to ball arithmetic: precision limits and unbounded exponents, aliased arguments and algebraic simplification, cached constants, and the use of hardware floating point in rigorous error bounds.
order: 4
---

The ball arithmetic modules follow the conventions of all FLINT modules
for types, naming, aliasing and memory management
([](guide/conventions)), and the general remarks on correctness, input
validation and size limits apply to them as well
([](guide/correctness)). The sections below describe the additional
conventions that are specific to ball arithmetic and easy to overlook.

## Precision limits {#limits}

Precisions are `slong` values and, like other sizes in FLINT
([](guide/correctness/limits)), are not checked against overflow. There
is no hard upper limit on the precision, but $2^{36} \approx 10^{11}$
bits on 64-bit systems and $2^{24}$ bits on 32-bit systems can be
considered safe. The limit on 64-bit systems comes from GMP, whose
integers, used internally by some algorithms, have at most about $2^{37}$
bits. Beware in particular of loops that keep doubling the precision
until a result is accurate enough: they need an upper bound, since the
result may never become accurate (see
[](guide/ball-arithmetic/using#guess-and-verify)). The minimum precision
is 2 bits.

These limits do not apply to exponents. The midpoints (`arf_t`) and
radii (`mag_t`) of balls have arbitrary-precision exponents, so
$2^{2^{100}}$ is an ordinary ball, and they do not apply to integers used
as scalars (as in `arb_mul_si`). They do apply to conversions whose
exact result has a size that depends on the exponent: converting
$2^{2^{100}}$ to an integer with `arf_get_fmpz` can do anything from
returning a wrong value to thrashing and crashing, and it is the
caller's responsibility not to try.

## Aliasing {#aliasing}

As everywhere in FLINT, inputs and outputs may be aliased unless
documented otherwise: in `arb_add(z, x, y, prec)` any two or all three
of `x`, `y` and `z` may be the same object. The rule concerns whole
objects of the same type ([](guide/conventions/types-and-naming#aliasing)),
and it matters in ball arithmetic because balls are often parts of
other objects. In `arb_poly_evaluate`, which computes $y = f(x)$, the
output `y` must not be a coefficient of `f`, although `x` and `y` may be
the same, and `x` may be a coefficient of `f`. The same holds for
`_arb_poly_evaluate`: an `arb_srcptr` array and an `arb_t` are not the
same type for aliasing purposes, even though a pointer to an array of
length 1 could otherwise be used as an `arb_t`.

## Aliased inputs and algebraic simplification {#simplification}

Identical input pointers give the implementation *permission* to
simplify algebraically. For example, `acb_mul(z, x, x, prec)` computes an
enclosure of

$$
\{ t^2 : t \in x \}
$$

using a squaring formula, and not of the generally larger set $\{t u : t
\in x, u \in x\}$. Consequently, if two values $x$ and $y$ are known to lie
in $[-1, 1]$ and an enclosure of $f(x, y)$ is wanted, it is a mistake to
create one ball for $[-1, 1]$ and call $f(x, x)$: two distinct balls,
each containing $[-1, 1]$, must be passed.

Simplification is allowed but not guaranteed. `arb_add(z, x, x, prec)`
and `arb_sub(z, x, x, prec)` treat their inputs as independent, and so
currently does `arb_mul(z, x, x, prec)` (and `arb_sqr`). To double a ball
use `arb_mul_2exp_si(z, x, 1)`, and to get zero use `arb_zero`.

## Cached constants {#caches}

Constants such as $\pi$, $\log 2$ and Euler's constant are read from
precomputed tables at low precision (up to 4608 bits for $\pi$) and
computed and cached at higher precision, so that later calls at the
same or a lower precision only round the cached value. The caches of
`arb_const_pi` and the other constants are per thread, and so are the
caches of Bernoulli numbers and of the logarithms of small primes used
by the elementary functions; some tables of [](mp_real) are shared
between threads.
None of this affects results, which are always rounded to the requested
precision. The caches are freed by `flint_cleanup` and
`flint_cleanup_master`; see [](guide/conventions/memory#cleanup) and,
for programs with several threads, [](guide/conventions/threading).

## Hardware floating point {#doubles}

The ball arithmetic modules use `double` arithmetic in two ways.

First, `double` arithmetic and `libm` functions such as `exp` and `log`
are used to choose parameters heuristically, for example the number of
terms of a series. This does not affect correctness: whenever an error
bound depends on such a parameter, the bound is computed rigorously. At
worst, a flawed floating-point implementation leads to inefficient
parameters.

Second, `double` arithmetic is used in some rigorous error bound
calculations, for example in `mag_t` arithmetic and in
`arb_sin_cos_wide`. This relies on the following assumptions, which hold
on all common platforms:

- `double` is the IEEE 754 binary64 format, with a 53-bit significand,
  an 11-bit exponent and encodings of infinities and NaN.
- The compiler does not perform unsafe floating-point optimisations such
  as reordering operations. These are off by default in GCC and Clang
  (and must not be enabled with options such as `-ffast-math`); some
  other compilers, such as the classic Intel compiler, enable them by
  default and must be told not to.
- Operations are not assumed to be correctly rounded (the x87 unit is not)
  or to round in a particular direction (the rounding mode may have been
  changed), only to have an error of at most 1.1 units in the last place.
- Underflow and overflow are not assumed to behave in any particular way:
  only `double` values in the normal range, and explicit infinities, are
  used.
- Transcendental `libm` functions are not used for rigorous bounds,
  since their error is not specified. Functions such as `ldexp` and
  `sqrt` are assumed to be correctly implemented. Where a bound needs a
  logarithm, the module uses its own `mag_d_log_upper_bound` and
  `mag_d_log_lower_bound`, which have proven error bounds.

Other parts of FLINT make stronger assumptions about `double`
arithmetic, for example the floating-point FFTs of [](fft_small) and the
`_precomp` functions of [](ulong_extras), so FLINT as a whole must not
be compiled with unsafe floating-point optimisations either.
