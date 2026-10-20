---
title: fmpz.h
header: fmpz.h
description: Arbitrary-precision integers with an inline representation of small values: arithmetic, division with all rounding conventions, GCDs, roots, modular arithmetic, Chinese remaindering, bit operations, primality testing and combinatorial functions.
category: Integers
order: 1
---

This module implements the integers $\Z$. An `fmpz_t` holds an integer of
any size and is the coefficient type of most of FLINT: polynomials,
matrices, rationals and number fields over $\Z$ are arrays of `fmpz`. A
value whose absolute value is at most `COEFF_MAX` $= 2^{\bits - 2} - 1$ is
stored directly in one machine word; larger values are stored in a GMP
`mpz` structure that the word points to. The switch between the two forms
is automatic, so `fmpz_t` behaves like GMP's `mpz_t`, but small values
cost no memory allocation and most operations on them run without a
function call into GMP. Large values use FLINT's own multiplication,
division and square root code where it is faster than GMP's.

Start with [](fmpz/getting-started), and read [](fmpz/representation)
before writing code that accesses the representation or creates large
arrays of integers. The naming, memory and aliasing conventions common to all FLINT modules
are described in [](guide/conventions/types-and-naming).

::: summary algorithms/integers/multiplication :::

## Types and conventions {#conventions}

::: type
fmpz
fmpz_t

An `fmpz` is an `slong`. If its two most significant bits are equal
($00$ or $11$), it is an integer $x$ with $|x| \le$ `COEFF_MAX`, stored
directly. If they are $01$, the remaining bits hold a pointer to an
`mpz_struct` shifted right by two bits, and the value is that of the
`mpz`. An `fmpz_t` is an array of one `fmpz`, so that it is passed by
reference like an `mpz_t`. A pointer `fmpz *` to an entry of an array can
be passed wherever an `fmpz_t` is expected.

::: details
Every `fmpz` must be initialised (`fmpz_init`, or a zero word) before use
and cleared afterwards (`fmpz_clear`), since a large value owns memory.
An array of integers is an array of `fmpz`, usually created with
`_fmpz_vec_init` and freed with `_fmpz_vec_clear`; see
[](fmpz/getting-started#arrays). The representation is described in
[](fmpz/representation).
:::
:::

::: macro
COEFF_MAX
COEFF_MIN

The largest and the smallest value stored inline: $2^{\bits - 2} - 1$ and
$-(2^{\bits - 2} - 1)$, that is, $\pm(2^{62} - 1)$ on 64-bit and
$\pm(2^{30} - 1)$ on 32-bit machines. The macro
`SMALL_FMPZ_BITCOUNT_MAX` is $\bits - 2$.
:::

::: macro
COEFF_IS_MPZ(f)

Nonzero if the `fmpz` value `f` (a word, not an `fmpz_t`) is a pointer to
an `mpz`, zero if it is a small value.
:::

::: function
mpz_ptr COEFF_TO_PTR(fmpz f)
fmpz PTR_TO_COEFF(mpz_ptr ptr)

Macros converting between an `fmpz` that holds a pointer and the
`mpz_ptr` itself. `PTR_TO_COEFF` shifts the pointer right by two bits and
sets the tag bits $01$; this relies on `mpz_struct`s being aligned to at
least four bytes. `COEFF_TO_PTR` must only be applied to values for which
`COEFF_IS_MPZ` is true.

::: details
`FMPZ_TO_ZZ` performs the same conversion for code that does not include
`gmp.h`.
:::
:::

::: type
zz_struct
zz_ptr
zz_srcptr

A structure with the same layout as GMP's `mpz_struct`, declared in
`fmpz_types.h`: the fields are `alloc` (the number of allocated limbs),
`size` (the number of limbs in use, negated for a negative value) and
`ptr` (the limbs, least significant first). It lets FLINT's headers
access the limbs of a large `fmpz` without including `gmp.h`.
:::

::: macro
FMPZ_TO_ZZ(x)

Converts an `fmpz` value `x` that holds a pointer (`COEFF_IS_MPZ(x)`
true) to a `zz_ptr`, like `COEFF_TO_PTR`.
:::

::: macro
MPZ_MIN_ALLOC

The minimum number of limbs allocated for an `mpz` owned by an `fmpz`,
currently 2. Every `mpz` handed out by `_fmpz_new_mpz` has at least this
many limbs, and every `mpz` given back to `_fmpz_clear_mpz` must still
have at least this many.
:::

The following rules apply to all functions of the module.

- *Aliasing.* Every output may be the same object as any input, unless
  stated otherwise. When a function has several outputs, such as a
  quotient and a remainder, the outputs must be distinct objects.
- *Canonical form.* A value that fits inline is always stored inline. Every
  function expects its inputs in this form and produces its outputs in
  it; only the internal functions documented at the end may leave a
  small value in an `mpz`. See [](fmpz/representation#canonical).
- *Errors.* Division by zero and other invalid arguments raise an
  exception through `flint_throw` (by default, a message and `abort`),
  unless stated otherwise. A few functions document preconditions that
  are not checked.
- *Functions with a leading underscore* are internal or have additional
  preconditions, stated in each case.

## Memory management {#memory}

::: function
void fmpz_init(fmpz_t f)

Initialises $f$ to zero. This only stores a zero word and allocates
nothing.
:::

::: function
void fmpz_init2(fmpz_t f, ulong limbs)

Initialises $f$ to zero with space for `limbs` limbs. If `limbs` is zero,
this is `fmpz_init`; otherwise an `mpz` is attached and grown to `limbs`
limbs, so that $f$ is temporarily not in canonical form (its value is 0 in
an `mpz`) until the first assignment. Only useful to avoid reallocations
when the final size is known.
:::

::: function
void fmpz_clear(fmpz_t f)

Releases the memory used by $f$. A large value gives its `mpz` back to the
memory manager (see [](fmpz/representation#memory)); a small value
needs no work.
:::

::: function
void fmpz_init_set(fmpz_t f, const fmpz_t g)
void fmpz_init_set_ui(fmpz_t f, ulong g)
void fmpz_init_set_si(fmpz_t f, slong g)

Initialise $f$ and set it to $g$.
:::

## Assignment and conversion {#conversion}

::: function
void fmpz_set(fmpz_t f, const fmpz_t g)

Sets $f$ to $g$.
:::

::: function
void fmpz_swap(fmpz_t f, fmpz_t g)

Swaps $f$ and $g$ by exchanging the two words; no limbs are copied.
:::

::: function
void fmpz_zero(fmpz_t f)
void fmpz_one(fmpz_t f)

Set $f$ to 0 and 1. A large value of $f$ is released.
:::

::: function
void fmpz_set_si(fmpz_t f, slong val)
void fmpz_set_ui(fmpz_t f, ulong val)
void fmpz_neg_ui(fmpz_t f, ulong val)

Set $f$ to `val`, and `fmpz_neg_ui` to $-$`val`.
:::

::: function
void fmpz_set_uiui(fmpz_t f, ulong hi, ulong lo)
void fmpz_neg_uiui(fmpz_t f, ulong hi, ulong lo)

Set $f$ to $h 2^{\bits} + l$, respectively to its negative, where $h$ and
$l$ are the unsigned words `hi` and `lo`.
:::

::: function
void fmpz_set_signed_uiui(fmpz_t r, ulong hi, ulong lo)
void fmpz_set_signed_uiuiui(fmpz_t r, ulong hi, ulong mid, ulong lo)

Set $r$ to the integer whose two's complement representation with
$2\bits$ (respectively $3\bits$) bits is `hi`, `lo` (respectively `hi`,
`mid`, `lo`), the most significant word first. These convert the results
of double-word signed arithmetic such as `smul_ppmm`.
:::

::: function
void fmpz_set_ui_array(fmpz_t out, const ulong * in, slong n)
void fmpz_set_signed_ui_array(fmpz_t out, const ulong * in, slong n)

Set `out` to the integer with the $n \ge 1$ limbs `in[0], ..., in[n-1]`,
least significant first: interpreted as an unsigned integer, or as a
two's complement integer of $n \bits$ bits.
:::

::: function
void fmpz_get_ui_array(ulong * out, slong n, const fmpz_t in)
void fmpz_get_signed_ui_array(ulong * out, slong n, const fmpz_t in)

Write the $n \ge 1$ limbs of `in` to `out`, least significant first. The
unsigned version requires $0 \le$ `in` $< 2^{n \bits}$; the signed version
writes `in` modulo $2^{n \bits}$, which is its two's complement
representation if it fits in $n \bits$ bits.
:::

::: function
void fmpz_set_mpn_large(fmpz_t z, nn_srcptr src, slong n, int negative)

Sets $z$ to the integer with the $n$ limbs `src`, negated if `negative` is
1. Requires $n \ge 2$ and `src[n-1]` nonzero, so that the value does not
fit inline; use `fmpz_set_ui` and `fmpz_neg_ui` for one limb, or
`fmpz_set_ui_array` in general.
:::

::: function
slong fmpz_get_si(const fmpz_t f)
ulong fmpz_get_ui(const fmpz_t f)

Return $f$ as an `slong`, respectively an `ulong`. The result is
unspecified if $f$ does not fit (use `fmpz_fits_si` and
`fmpz_abs_fits_ui` to check); currently `fmpz_get_ui` returns the low limb
of $|f|$, so a negative $f$ is not reduced modulo $2^{\bits}$.
:::

::: function
void fmpz_get_uiui(ulong * hi, ulong * low, const fmpz_t f)

Sets `*hi` and `*low` to the high and low limbs of $f$, for
$0 \le f < 2^{2\bits}$.
:::

::: function
void fmpz_get_signed_uiui(ulong * hi, ulong * lo, const fmpz_t x)

Sets `*hi` and `*lo` to the two limbs of $x$ modulo $2^{2\bits}$, its two's
complement representation if $x$ fits in $2\bits$ bits.
:::

::: function
ulong fmpz_get_nmod(const fmpz_t f, nmod_t mod)

Returns $f \bmod n$ in $[0, n)$ for the modulus $n$ of `mod`, for any $f$.
:::

::: function
double fmpz_get_d(const fmpz_t f)

Returns $f$ as a `double`, rounded towards zero if it is not exactly
representable. The result is undefined if $|f|$ exceeds the range of
`double`; use `fmpz_get_d_2exp` for large values.
:::

::: function
void fmpz_set_d(fmpz_t f, double c)

Sets $f$ to $c$ rounded towards zero. The result is undefined if $c$ is
infinite or a NaN.
:::

::: function
double fmpz_get_d_2exp(slong * exp, const fmpz_t f)

Returns a `double` $d$ with $0.5 \le |d| < 1$ and sets `*exp` to $e$ such
that $d 2^e$ approximates $f$; $d$ is $f 2^{-e}$ truncated to 53 bits. For
$f = 0$ it returns 0 and sets `*exp` to 0.
:::

::: function
void fmpz_set_d_2exp(fmpz_t f, double d, slong exp)

Sets $f$ to $d 2^{\mathit{exp}}$ rounded towards zero, so that
$0.75 \cdot 2^0$ gives 0.
:::

::: function
void fmpz_get_mpz(mpz_t x, const fmpz_t f)
void fmpz_set_mpz(fmpz_t f, const mpz_t x)

Convert between `fmpz_t` and an initialised GMP `mpz_t`, copying the
value. Declared only if `gmp.h` is included before the FLINT headers.
:::

::: function
void fmpz_get_mpf(mpf_t x, const fmpz_t f)
void fmpz_set_mpf(fmpz_t f, const mpf_t x)

Convert between `fmpz_t` and GMP's `mpf_t`; `fmpz_set_mpf` truncates a
fractional value towards zero. Declared only if `gmp.h` is included first.
:::

::: function
void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)

Sets the MPFR number $x$ to $f$, rounded in the direction `rnd` to the
precision of $x$. Declared only if `mpfr.h` is included before the FLINT
headers.
:::

::: function
int fmpz_get_mpn(nn_ptr * n, fmpz_t n_in)

Allocates an array with `flint_malloc`, stores the limbs of $|n_{in}|$ in
it, sets `*n` to it and returns the number of limbs (at least 1). The
caller frees the array with `flint_free`. For $n_{in} = 0$ the array has
one limb, which is 0.
:::

## Converting to and from GMP without copying {#readonly}

These functions give GMP read access to an `fmpz_t` and FLINT read
access to an `mpz_t`, without copying limbs. The guide
[](fmpz/representation#gmp) shows how to use them. All require `gmp.h` to
be included before the FLINT headers.

::: function
void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f)
void flint_mpz_clear_readonly(mpz_t z)

`flint_mpz_init_set_readonly` sets the uninitialised `mpz_t` $z$ to the
value of $f$, sharing the limbs of $f$ if it is large; $z$ must not be
modified, $f$ must not change while $z$ is in use, and $z$ is released with
`flint_mpz_clear_readonly`.
:::

::: function
void fmpz_init_set_readonly(fmpz_t f, const mpz_t z)
void fmpz_clear_readonly(fmpz_t f)

`fmpz_init_set_readonly` sets the uninitialised `fmpz_t` $f$ to the value
of $z$, sharing the limbs of $z$ if the value is large; $f$ must not be
modified, $z$ must not change while $f$ is in use, and $f$ is released
with `fmpz_clear_readonly`.
:::

## Input and output {#io}

The function page [](fmpz/strings) gives the details and an example.

::: function
char * fmpz_get_str(char * str, int b, const fmpz_t f)

Returns the representation of $f$ in base $b$, $2 \le b \le 62$, with a
leading minus sign if $f < 0$. If `str` is `NULL`, the string is allocated
with `flint_malloc` and must be freed with `flint_free`; otherwise it is
written to `str`, which must have room for
`fmpz_sizeinbase(f, b) + 2` characters.
:::

::: function
int fmpz_set_str(fmpz_t f, const char * str, int b)

Sets $f$ to the integer written in base $b$, $2 \le b \le 62$, in the
null-terminated string `str`, and returns 0. White space before and
after the number and a minus sign are accepted. Returns $-1$ if the string
is not a valid integer, in which case $f$ is undefined.
:::

::: function
int fmpz_fprint(FILE * file, const fmpz_t x)
int fmpz_print(const fmpz_t x)

Print $x$ in decimal to `file` or to `stdout`, with a minus sign if
negative and no newline. Return the number of characters written, or a
nonpositive value on error.
:::

::: function
int fmpz_fread(FILE * file, fmpz_t f)
int fmpz_read(fmpz_t f)

Read a decimal integer (optional white space and minus sign, then digits)
from `file` or from `stdin` into $f$. Return a positive value on success
and 0 on failure, following `mpz_inp_str`.
:::

::: function
size_t fmpz_out_raw(FILE * fout, const fmpz_t x)
size_t fmpz_inp_raw(fmpz_t x, FILE * fin)

Write $x$ to, respectively read it from, a stream in GMP's portable binary
format: four bytes holding the signed number of bytes, then the bytes of
$|x|$, most significant first. Files are interchangeable with
`mpz_out_raw` and `mpz_inp_raw`. Return the number of bytes written or
read, and 0 on error.
:::

## Random integers {#random}

All functions take a `flint_rand_t`; see [](flint). They are meant for
testing and are not suitable for cryptography.

::: function
void fmpz_randbits_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
void fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)

Set $f$ to a uniformly random integer with exactly `bits` bits (its top
bit set), nonnegative for the unsigned version and with a random sign
otherwise. For `bits` $= 0$ the result is 0.
:::

::: function
void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
void fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)

Set $f$ to a random integer whose absolute value has a random number of
bits between 0 and `bits`, nonnegative for the unsigned version. The
values are biased towards sparse binary representations, 0, 1 and
`COEFF_MAX`, which exercise special cases in tests.
:::

::: function
void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)

As `fmpz_randtest`, but never 0. Requires `bits` $\ge 1$.
:::

::: function
void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m)

Sets $f$ to a uniformly random integer in $[0, m)$ for $m > 0$, and in
$(m, 0]$ for $m < 0$.
:::

::: function
void fmpz_randm_not_zero(fmpz_t f, flint_rand_t state, const fmpz_t m)

Sets $f$ to a uniformly random integer in $[1, m)$. Requires $m \ge 2$.
:::

::: function
void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m)
void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m)

Set $f$ to a random integer in $[0, m)$, respectively $(-m/2, m/2]$, for
$m > 0$, with an increased probability of values near the ends of the
range (and near zero for the signed version).
:::

::: function
void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved)

Sets $f$ to a random prime with exactly `bits` $\ge 2$ bits: the next
prime after a random integer of that size, retried until it has the
right size, so the distribution is not exactly uniform over primes. If
`proved` is nonzero the result is proved prime, otherwise it is a
probable prime as in `fmpz_nextprime`.
:::

## Properties and comparison {#properties}

::: function
int fmpz_sgn(const fmpz_t f)

Returns $-1$, $0$ or $1$ according to the sign of $f$.
:::

::: function
int fmpz_is_zero(const fmpz_t f)
int fmpz_is_one(const fmpz_t f)
int fmpz_is_pm1(const fmpz_t f)

Return whether $f = 0$, $f = 1$ and $f = \pm 1$.
:::

::: function
int fmpz_is_even(const fmpz_t f)
int fmpz_is_odd(const fmpz_t f)

Return whether $f$ is even, respectively odd.
:::

::: function
int fmpz_equal(const fmpz_t f, const fmpz_t g)
int fmpz_equal_si(const fmpz_t f, slong g)
int fmpz_equal_ui(const fmpz_t f, ulong g)

Return 1 if $f = g$ and 0 otherwise.
:::

::: function
int fmpz_cmp(const fmpz_t f, const fmpz_t g)
int fmpz_cmp_ui(const fmpz_t f, ulong g)
int fmpz_cmp_si(const fmpz_t f, slong g)

Return a negative value, zero or a positive value according to whether
$f < g$, $f = g$ or $f > g$.
:::

::: function
int fmpz_cmpabs(const fmpz_t f, const fmpz_t g)
int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g)

Compare $|f|$ with $|g|$, respectively with $|2g|$, returning a negative
value, zero or a positive value.
:::

::: function
flint_bitcnt_t fmpz_bits(const fmpz_t f)
slong fmpz_size(const fmpz_t f)

Return the number of bits, respectively limbs, of $|f|$; both are 0 for
$f = 0$.
:::

::: function
size_t fmpz_sizeinbase(const fmpz_t f, int b)

Returns the number of digits of $|f|$ in base $b$, $2 \le b \le 62$, and
1 for $f = 0$. Unlike GMP's `mpz_sizeinbase`, the result is exact for all
bases.
:::

::: function
flint_bitcnt_t fmpz_val2(const fmpz_t x)

Returns the exponent of the largest power of 2 dividing $x$ (the number of
trailing zero bits), and 0 for $x = 0$.
:::

::: function
int fmpz_abs_fits_ui(const fmpz_t f)
int fmpz_fits_si(const fmpz_t f)

Return whether $|f|$ fits in an `ulong`, respectively whether $f$ fits in
an `slong`.
:::

::: function
ulong fmpz_abs_lbound_ui_2exp(slong * exp, const fmpz_t x, int bits)
ulong fmpz_abs_ubound_ui_2exp(slong * exp, const fmpz_t x, int bits)

For $x \ne 0$ and $1 \le$ `bits` $\le \bits$, return a mantissa $m$ with
exactly `bits` bits and set `*exp` to $e$ such that $m 2^e \le |x|$
(lower bound, $m$ is $|x|$ truncated), respectively $m 2^e \ge |x|$
(upper bound, $m$ is $|x|$ rounded up, which may increase $e$ by one if
the rounding carries).
:::

::: function
double fmpz_dlog(const fmpz_t x)

Returns an approximation of $\log x$ for $x > 0$, with a relative error of
a few units in the last place of a `double`.
:::

::: function
slong fmpz_flog(const fmpz_t x, const fmpz_t b)
slong fmpz_flog_ui(const fmpz_t x, ulong b)
slong fmpz_clog(const fmpz_t x, const fmpz_t b)
slong fmpz_clog_ui(const fmpz_t x, ulong b)

Return $\lfloor \log_b x \rfloor$ and $\lceil \log_b x \rceil$, computed
exactly. Require $x \ge 1$ and $b \ge 2$.
:::

## Addition and multiplication {#arithmetic}

::: function
void fmpz_neg(fmpz_t f1, const fmpz_t f2)
void fmpz_abs(fmpz_t f1, const fmpz_t f2)

Set $f_1$ to $-f_2$, respectively $|f_2|$.
:::

::: function
void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x)
void fmpz_add_si(fmpz_t f, const fmpz_t g, slong x)
void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x)
void fmpz_sub_si(fmpz_t f, const fmpz_t g, slong x)

Set $f$ to $g + h$, respectively $g - h$ (with $x$ in place of $h$).
:::

::: function
void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x)
void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x)

Set $f$ to $g h$ (with $x$ in place of $h$). Two small operands are
multiplied inline; large ones go to `flint_mpn_mul` or `flint_mpn_sqr`,
which switch from basecase code through Karatsuba and GMP's Toom–Cook
algorithms to FFT multiplication ([](algorithms/integers/multiplication)).
:::

::: function
void fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, ulong h1, ulong h2)

Sets $f$ to $g h_1 h_2$.
:::

::: function
void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x)
void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong x)
void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x)
void fmpz_submul_si(fmpz_t f, const fmpz_t g, slong x)

Set $f$ to $f + g h$, respectively $f - g h$ (with $x$ in place of $h$),
without a temporary for the product when the operands are small.
:::

::: function
void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)
void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)

Set $f$ to $ab + cd$, respectively $ab - cd$.
:::

::: function
void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_one_2exp(fmpz_t f, ulong exp)

Set $f$ to $g 2^{\mathit{exp}}$, respectively $2^{\mathit{exp}}$.
`fmpz_mul_2exp` assumes that $\mathit{exp} + \bits$ does not overflow an
`ulong`.
:::

::: function
void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp)
void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slong x, ulong exp)

Set $f$ to $g h / 2^{\mathit{exp}}$ (with $x$ in place of $h$), rounded
towards zero.
:::

## Powers and roots {#roots}

The function page [](fmpz/roots) describes the root functions and their
algorithms.

::: function
void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)

Set $f$ to $g^{\mathit{exp}}$, respectively $x$ to $b^e$, with
$0^0 = 1$.
:::

::: function
int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e)

Sets $a$ to $b^e$ with $0^0 = 1$ and returns 1. For $|b| \ge 2$ and an
exponent $e$ that does not fit in an `slong` (the result would not fit
in memory), returns 0 and leaves $a$ unchanged. Throws an exception if
$e < 0$.
:::

::: function
void fmpz_sqrt(fmpz_t f, const fmpz_t g)
void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g)

Set $f$ to $\lfloor \sqrt g \rfloor$ and $r$ to $g - f^2$, for $g \ge 0$.
Throw an exception if $g < 0$. In `fmpz_sqrtrem`, $f$ and $r$ must be
distinct objects.
:::

::: function
int fmpz_perfect_sqrt(fmpz_t s, const fmpz_t a)

If $a$ is a perfect square, sets $s$ to $\sqrt a \ge 0$ and returns 1;
otherwise returns 0 and leaves $s$ undefined. Faster than `fmpz_is_square`
followed by `fmpz_sqrt`.
:::

::: function
int fmpz_is_square(const fmpz_t f)

Returns whether $f$ is a perfect square ($f \ge 0$).
:::

::: function
int fmpz_root(fmpz_t r, const fmpz_t f, slong n)

Sets $r$ to the $n$-th root of $f$ rounded towards zero and returns 1 if
it is exact and 0 otherwise. Requires $n \ge 1$, and $f \ge 0$ if $n$ is
even; otherwise an exception is raised.
:::

::: function
int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f)

If $f = r^k$ for some $k \ge 2$, sets `root` to $r$ and returns $k$;
otherwise returns 0. Negative $f$ are allowed (then $k$ is odd). The
exponent is not necessarily the largest one: $2^{120}$ currently gives
$k = 2$, and $-2^{120}$ gives $k = 3$. By convention 0 and 1 give $k = 2$
and $-1$ gives $k = 1$.
:::

## Division {#division}

Division comes in four rounding conventions: `tdiv` rounds the quotient
towards zero, `fdiv` towards $-\infty$, `cdiv` towards $+\infty$ and
`ndiv` to the nearest integer. The function page [](fmpz/division)
explains them, and [](algorithms/integers/division) the algorithms.

::: function
void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)

Set $f$ to the quotient $q$ of $g$ by $h$, rounded towards zero, $-\infty$
and $+\infty$ respectively, and $s$ to the remainder $g - qh$. The
remainder has the sign of $g$ (`tdiv`), of $h$ (`fdiv`) or of $-h$
(`cdiv`), and $|s| < |h|$. Throw an exception if $h = 0$.
:::

::: function
void fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)

Sets $q$ to $a/b$ rounded to the nearest integer, ties to even, and $r$ to
$a - qb$, so that $|r| \le |b|/2$. Throws an exception if $b = 0$.
:::

::: function
void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h)

Set $f$ to the quotient of $g$ by $h$, rounded towards zero, $-\infty$
and $+\infty$ respectively. Throw an exception if $h = 0$.
:::

::: function
void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)

Sets $f$ to the remainder of `fmpz_fdiv_qr`, which has the sign of $h$.
Throws an exception if $h = 0$.
:::

::: function
void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_tdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
void fmpz_cdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)

Set $f$ to the quotient, respectively the remainder, of $g$ by
$2^{\mathit{exp}}$ with the rounding of the prefix. `fmpz_fdiv_q_2exp` is
an arithmetic right shift and `fmpz_fdiv_r_2exp` the reduction to
$[0, 2^{\mathit{exp}})$.
:::

::: function
ulong fmpz_tdiv_ui(const fmpz_t g, ulong h)
ulong fmpz_fdiv_ui(const fmpz_t g, ulong h)
ulong fmpz_cdiv_ui(const fmpz_t g, ulong h)

Return the absolute value of the remainder of $g$ by $h$ with the
rounding of the prefix. `fmpz_fdiv_ui` thus returns $g \bmod h$ in
$[0, h)$. Throw an exception if $h = 0$.
:::

::: function
void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h)
void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h)

Set $f$ to $g / h$, assuming that $h$ divides $g$. This is not checked;
if the division is not exact the result is meaningless. Faster than the
other divisions, since it uses 2-adic (Hensel) division. Throw an
exception if $h = 0$.
:::

::: function
void fmpz_divexact2_uiui(fmpz_t f, const fmpz_t g, ulong h1, ulong h2)

Sets $f$ to $g / (h_1 h_2)$, assuming that $h_1 h_2$ divides $g$. Throws
an exception if $h_1$ or $h_2$ is 0.
:::

::: function
int fmpz_div(fmpz_t q, const fmpz_t g, const fmpz_t h)
int fmpz_div_ui(fmpz_t q, const fmpz_t g, ulong h)
int fmpz_div_si(fmpz_t q, const fmpz_t g, slong h)

Checked exact division: if $h$ divides $g$, set $q = g/h$ and return 1;
otherwise set $q = 0$ and return 0. Division by zero counts as exact only
for $g = 0$ (with $q = 0$).
:::

::: function
int fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h)

Older name for `fmpz_div`.
:::

::: function
int fmpz_divisible(const fmpz_t f, const fmpz_t g)

Returns whether $g$ divides $f$, that is, whether $f = qg$ for an integer
$q$. Every $g$ divides 0, and 0 divides only 0.
:::

::: function
int fmpz_divisible_ui(const fmpz_t f, ulong g)
int fmpz_divisible_si(const fmpz_t f, slong g)

Return whether $g$ divides $f$. As for `fmpz_divisible`, every $g$
divides 0, and 0 divides only 0.
:::

::: type
fmpz_preinvn_struct
fmpz_preinvn_t

A precomputed inverse of a divisor $h$ for `fmpz_fdiv_qr_preinvn`: the
fields are `dinv` (an $n$-limb approximation of the reciprocal of $h'$,
computed by `flint_mpn_preinvn`, where $h'$ is $|h|$ shifted left by
`norm` bits to have its top bit set) and `n` (the number of limbs of
$h$).
:::

::: function
void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f)
void fmpz_preinvn_clear(fmpz_preinvn_t inv)

Compute a precomputed inverse of $f \ne 0$, respectively release it.
:::

::: function
void fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t inv)
void fmpz_fdiv_r_preinvn(fmpz_t f, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t inv)

As `fmpz_fdiv_qr` and `fmpz_fdiv_r`, using the precomputed inverse `inv`
of $h$, for many divisions by the same multi-limb $h$.
:::

## Modular arithmetic {#modular}

The function page [](fmpz/modular) describes inversion and
exponentiation.

::: function
void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)

Sets $f$ to $g \bmod h$ in $[0, |h|)$. Requires $h \ne 0$, which is not
checked (a zero divisor crashes with a division by zero).
:::

::: function
ulong fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h)

Sets $f$ to $g \bmod h$ in $[0, h)$ and also returns it. Throws an
exception if $h = 0$.
:::

::: function
void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h)

Sets $f$ to the symmetric remainder $y \equiv g \pmod h$ with
$-|h|/2 < y \le |h|/2$. Requires $h \ne 0$ (not checked).
:::

::: function
void fmpz_set_ui_smod(fmpz_t f, ulong x, ulong m)

Sets $f$ to the symmetric remainder $y \equiv x \pmod m$ with
$-m/2 < y \le m/2$, for $0 \le x < m$ (not checked).
:::

::: function
void fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)

Sets $r$ to $-a$ modulo `mod` for $0 \le a <$ `mod` (not checked).
:::

::: function
int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h)

If $\gcd(g, h) = 1$, sets $f$ to the inverse of $g$ modulo $h$ in
$[0, |h|)$ and returns 1; otherwise returns 0 and $f$ is undefined. Every
$g$ is invertible modulo $h = \pm 1$, with inverse 0. Throws an exception
if $h = 0$.
:::

::: function
void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong exp, const fmpz_t m)
void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)

Set $f$ to $g^e \bmod m$ in $[0, m)$, for $m \ge 1$ (an exception is
raised otherwise). For $m = 1$ the result is 0, and otherwise $g^0 = 1$.
`fmpz_powm` accepts $e < 0$ if $g$ is invertible modulo $m$ and raises an
exception if it is not.
:::

::: function
int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p)

For a prime $p$: if $a$ is a square modulo $p$, sets $b$ to a square root
and returns 1, otherwise sets $b$ to 0 and returns 0. Primality of $p$ is
not checked; for composite $p$ a return value of 1 means nothing.
:::

::: function
slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)

Sets `rop` to `op` with all factors $f$ removed and returns their number,
the largest $v$ with $f^v \mid$ `op`. Requires $f \ge 2$ (an exception is
raised otherwise). For `op` $= 0$, `rop` is set to 0 and 0 is returned.
:::

::: function
int fmpz_jacobi(const fmpz_t a, const fmpz_t p)

Returns the Jacobi symbol $\left(\frac{a}{p}\right)$ for any $a$ and odd
$p > 0$. The conditions on $p$ are not checked.
:::

::: function
int fmpz_kronecker(const fmpz_t a, const fmpz_t n)

Returns the Kronecker symbol $\left(\frac{a}{n}\right)$ for any $a$ and
$n$.
:::

::: function
void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n)

Sets the outputs so that the solutions $x \in [0, n)$ of $ax \equiv b
\pmod n$ are exactly $\mathit{xstart} + i \cdot \mathit{xstride}$ for $0
\le i < \mathit{xlength}$; there are none if `xlength` is 0. Requires
$n > 0$; the outputs must be distinct objects.
:::

### Arithmetic modulo powers of two {#mod-2exp}

These functions use 2-adic Newton iteration on the limbs
([](algorithms/integers/division#hensel)). Results are reduced to
$[0, 2^N)$.

::: function
void fmpz_invmod_2exp(fmpz_t rop, const fmpz_t op, flint_bitcnt_t N)

Sets `rop` to the inverse of `op` modulo $2^N$. Throws an exception if
`op` is even.
:::

::: function
void fmpz_divmod_2exp(fmpz_t rop, const fmpz_t a, const fmpz_t b, flint_bitcnt_t N)

Hensel division: sets `rop` to the unique $q \in [0, 2^N)$ with
$qb \equiv a \pmod{2^N}$. Throws an exception if $b$ is even.
:::

::: function
int fmpz_sqrtmod_2exp(fmpz_t rop, const fmpz_t op, flint_bitcnt_t N)
int fmpz_rsqrtmod_2exp(fmpz_t rop, const fmpz_t op, flint_bitcnt_t N)

If `op` is an odd square modulo $2^N$ (for $N \ge 3$: `op` $\equiv 1
\pmod 8$), set `rop` to a square root $y$, respectively an inverse square
root ($\mathit{op}\, y^2 \equiv 1$), modulo $2^N$, and return 1;
otherwise return 0. The root returned is the one congruent to 1 modulo
4, reduced modulo $2^N$. Even `op` is not supported (0 is returned).
:::

## Greatest common divisor {#gcd}

The function page [](fmpz/gcd) explains the extended GCD variants and
their normalisations; [](algorithms/integers/gcd) describes the
algorithms.

::: function
void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
void fmpz_gcd_ui(fmpz_t res, const fmpz_t a, ulong b)

Set $f$ to $\gcd(g, h) \ge 0$, with $\gcd(0, 0) = 0$.
:::

::: function
void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c)

Sets $f$ to $\gcd(a, b, c) \ge 0$. Faster than two calls to `fmpz_gcd`
when one operand is small or the GCD of the two shortest operands is 1.
:::

::: function
void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)

Sets $f$ to the least common multiple of $g$ and $h$, which is
nonnegative, and 0 if $g$ or $h$ is 0.
:::

::: function
void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)

Sets $d = \gcd(f, g) \ge 0$ and cofactors $a$, $b$ with $af + bg = d$.
The outputs must be distinct objects.
:::

::: function
void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)

As `fmpz_xgcd`, but the cofactors are the canonical Bézout coefficients
defined on [](fmpz/gcd): the minimal ones, with fixed choices in the
degenerate cases. Usually faster than `fmpz_xgcd`.
:::

::: function
void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)

For $0 \le f < g$, sets $d = \gcd(f, g)$ and $a$ to a cofactor in
$[0, g)$ with $af \equiv d \pmod g$, so that $a = f^{-1} \bmod g$ if
$d = 1$. For $f = 0$ it sets $d = g$ and $a = 0$. The outputs must be
distinct objects.
:::

::: function
void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, const fmpz_t L)

Runs the extended Euclidean algorithm on $r_2 > r_1 > 0$ until the
remainder drops to at most $L$, with Lehmer's acceleration. On return
$r_2$ and $r_1$ are the last two remainders and `co2`, `co1` the
corresponding cofactors, with `co2`$\cdot r_1 -$ `co1`$\cdot r_2 = \pm
r_{2,\mathrm{orig}}$. No aliasing is allowed.
:::

## Chinese remaindering {#crt}

The function page [](fmpz/crt) compares the CRT interfaces and has
examples; [](algorithms/integers/binary-splitting#remainder-trees)
describes the product and remainder trees used for many moduli.

::: function
void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2, ulong m2, int sign)
void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2, int sign)

Set `out` to the unique $x \equiv r_1 \pmod{m_1}$, $x \equiv r_2 \pmod
{m_2}$ in $[0, m_1 m_2)$ if `sign` is 0, and in $(-m_1 m_2/2, m_1
m_2/2]$ otherwise. Require coprime $m_1, m_2 > 0$, $-m_1 \le r_1 < m_1$
and $0 \le r_2 < m_2$; only the coprimality is checked (by `fmpz_CRT`,
which raises an exception).
:::

::: function
void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2, ulong m2, ulong m2inv, const fmpz_t m1m2, ulong c, int sign)

As `fmpz_CRT_ui` with precomputed data: `m2inv` $=$
`n_preinvert_limb(m2)`, `m1m2` $= m_1 m_2$ and $c = m_1^{-1} \bmod m_2$.
Used when the same pair of moduli recurs.
:::

::: type
fmpz_comb_struct
fmpz_comb_t
fmpz_comb_temp_struct
fmpz_comb_temp_t

Precomputed data for reduction modulo, and reconstruction from, a fixed
list of word-size moduli (`fmpz_comb_t`), and the scratch space for one
call (`fmpz_comb_temp_t`). An `fmpz_comb_t` is a thin wrapper around the
`flint_mpn_crt_t` of `mpn_extras`.
:::

::: macro
FMPZ_COMB_MOD
FMPZ_COMB_CRT

Flags for `fmpz_comb_init2` selecting the precomputation for
`fmpz_multi_mod_ui` and `fmpz_multi_CRT_ui` respectively.
:::

::: function
void fmpz_comb_init(fmpz_comb_t C, nn_srcptr primes, slong num_primes)
void fmpz_comb_init2(fmpz_comb_t C, nn_srcptr primes, slong num_primes, int flags)
void fmpz_comb_clear(fmpz_comb_t C)

Initialise `C` for the $\mathit{num\_primes} \ge 1$ pairwise coprime
nonzero word-size moduli in `primes` (they need not be prime; the array
is copied), respectively release it. `fmpz_comb_init2` precomputes only
what the operations in `flags` (`FMPZ_COMB_MOD`, `FMPZ_COMB_CRT` or both)
need.
:::

::: function
void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C)
void fmpz_comb_temp_clear(fmpz_comb_temp_t CT)

Allocate, respectively free, scratch space for calls with `C`. A scratch
space is not thread-safe; each thread needs its own.
:::

::: function
void fmpz_multi_mod_ui(ulong * out, const fmpz_t in, const fmpz_comb_t C, fmpz_comb_temp_t CT)

Sets `out[i]` to `in` modulo the $i$-th modulus of `C`, in $[0, m_i)$,
for any `in`.
:::

::: function
void fmpz_multi_CRT_ui(fmpz_t output, nn_srcptr residues, const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign)

Sets `output` to the integer congruent to `residues[i]` $\in [0, m_i)$
modulo each modulus of `comb`, in $[0, M)$ if `sign` is 0 and in
$(-M/2, M/2]$ otherwise, where $M$ is the product of the moduli.
:::

::: function
void fmpz_multi_mod_ui_once(ulong * out, const fmpz_t in, nn_srcptr primes, slong num_primes)
void fmpz_multi_CRT_ui_once(fmpz_t output, fmpz_t prod, nn_srcptr residues, nn_srcptr primes, slong num_primes, int sign)

One-shot versions of `fmpz_multi_mod_ui` and `fmpz_multi_CRT_ui` that do
not need an `fmpz_comb_t`. If `prod` is not `NULL`, the CRT version also
sets it to the product of the moduli.
:::

::: type
fmpz_multi_CRT_struct
fmpz_multi_CRT_t

Precomputed data for Chinese remaindering with a fixed list of moduli of
any size.
:::

::: function
void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT)
void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P)

Initialise `CRT`, respectively release `P`.
:::

::: function
int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT, const fmpz * moduli, slong len)

Prepares `CRT` for the $\mathit{len} \ge 1$ moduli (signs are ignored)
and returns 1 if they are usable: one nonzero modulus, or several that are
pairwise coprime and none of which is $0$ or $\pm 1$. Otherwise returns 0,
and later reconstructions give undefined results.
:::

::: function
void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign)

Sets `output` to the integer congruent to `inputs[i]` modulo the $i$-th
modulus of `P`, for any values `inputs[i]`: in $[0, M)$ if `sign` is 0
and in $(-M/2, M/2]$ otherwise, where $M$ is the product of the absolute
values of the moduli.
:::

::: function
void _fmpz_multi_CRT_precomp(fmpz * outputs, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign)

As `fmpz_multi_CRT_precomp`, with the result in `outputs[0]` and the
other entries of the array, which must have `P->localsize` initialised
entries, used as scratch space.
:::

::: function
int fmpz_multi_CRT(fmpz_t output, const fmpz * moduli, const fmpz * values, slong len, int sign)

Combines `fmpz_multi_CRT_precompute` and `fmpz_multi_CRT_precomp` and
returns the value of the former.
:::

::: type
fmpz_multi_mod_struct
fmpz_multi_mod_t

Precomputed data for reducing an integer modulo a fixed list of moduli of
any size.
:::

::: function
void fmpz_multi_mod_init(fmpz_multi_mod_t P)
void fmpz_multi_mod_clear(fmpz_multi_mod_t P)

Initialise, respectively release, `P`.
:::

::: function
int fmpz_multi_mod_precompute(fmpz_multi_mod_t P, const fmpz * f, slong r)

Prepares `P` for reduction modulo the $r \ge 1$ moduli `f` (signs are
ignored) and returns 1 if none of them is zero, 0 otherwise.
:::

::: function
void fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign)
void _fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign, fmpz * tmp)

Set `outputs[i]` to `input` modulo the $i$-th modulus $m_i$ of `P`: in
$[0, |m_i|)$ if `sign` is 0 and in $(-|m_i|/2, |m_i|/2]$ otherwise. The
underscore version takes scratch space `tmp` of `P->localsize`
initialised entries.
:::

## Bit operations {#bits}

Bit operations act on the infinite two's complement representation, as in
GMP: a negative integer has infinitely many leading one bits. The function
page [](fmpz/bits) has examples.

::: function
int fmpz_tstbit(const fmpz_t f, ulong i)

Returns bit $i$ of $f$ (0 or 1).
:::

::: function
void fmpz_setbit(fmpz_t f, ulong i)
void fmpz_clrbit(fmpz_t f, ulong i)
void fmpz_combit(fmpz_t f, ulong i)

Set, clear and complement bit $i$ of $f$.
:::

::: function
void fmpz_complement(fmpz_t r, const fmpz_t f)

Sets $r$ to the one's complement $-f - 1$ of $f$.
:::

::: function
void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b)
void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b)
void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b)

Set $r$ to the bitwise and, inclusive or and exclusive or of $a$ and $b$.
:::

::: function
flint_bitcnt_t fmpz_popcnt(const fmpz_t c)

Returns the number of one bits of $c \ge 0$. For negative $c$ (infinitely
many one bits) it returns 0.
:::

::: function
int fmpz_bit_pack(nn_ptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, const fmpz_t coeff, int negate, int borrow)
int fmpz_bit_unpack(fmpz_t coeff, nn_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, int negate, int borrow)
void fmpz_bit_unpack_unsigned(fmpz_t coeff, nn_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits)

Pack an integer into, respectively extract it from, a field of `bits`
bits starting at bit `shift` $<\bits$ of the limb array `arr`. Signed
fields are stored in two's complement with a borrow that propagates to
the next field. These are the building blocks of Kronecker substitution
in `fmpz_poly`.

::: details
`fmpz_bit_pack` adds `coeff` (negated first if `negate` is $-1$, minus
`borrow` $\in \{0, 1\}$) to the field, assuming that the bits of `arr`
from position `shift` on are zero. It returns 1 if the packed value is
negative, which is the borrow to pass when packing the next field.

`fmpz_bit_unpack` sets `coeff` to the field, read as a signed value of
`bits` bits, plus `borrow`, negated if `negate` is $-1$, and returns 1 if
the field was negative, which is the borrow to add to the next field.
`fmpz_bit_unpack_unsigned` reads the field as an unsigned value.
:::
:::

## Combinatorial functions {#combinatorics}

[](algorithms/integers/binary-splitting) explains the product trees behind
these functions.

::: function
void fmpz_fac_ui(fmpz_t f, ulong n)

Sets $f$ to $n!$, from a table for $n \le 20$ (on 64-bit machines) and
with GMP's `mpz_fac_ui` otherwise.
:::

::: function
void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k)

Sets `res` to the binomial coefficient $\binom{n}{k}$, which is 0 for
$k > n$.
:::

::: function
void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n)
void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n)

Set $r$ to the rising factorial $x(x+1)\cdots(x+n-1)$, with value 1 for
$n = 0$. Any sign of $x$ is allowed.
:::

::: function
void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b)

Sets $r$ to $(x+a)(x+a+1)\cdots(x+b-1)$ for $x \ge 0$ and $b > a$, by
binary splitting with several factors packed per word at the leaves.
:::

::: function
void fmpz_fib_ui(fmpz_t f, ulong n)

Sets $f$ to the Fibonacci number $F_n$.
:::

::: function
void fmpz_primorial(fmpz_t res, ulong n)

Sets `res` to the primorial $n\#$, the product of the primes $p \le n$,
by a balanced product tree over word-size factors.
:::

::: function
void fmpz_euler_phi(fmpz_t res, const fmpz_t n)

Sets `res` to Euler's totient $\varphi(n)$, the number of
$1 \le k \le n$ with $\gcd(k, n) = 1$, and to 0 for $n \le 0$. Word-size
$n$ use `n_euler_phi`; larger $n$ are factored. The function
`fmpz_factor_euler_phi` of `fmpz_factor` takes a factorisation.
:::

::: function
int fmpz_moebius_mu(const fmpz_t n)

Returns the Möbius function $\mu(|n|)$: 0 if $n$ has a square factor
$> 1$, otherwise $(-1)^k$ for $k$ prime factors; $\mu(0) = 0$. Factors
$n$ unless it fits in a word.
:::

::: function
void fmpz_divisor_sigma(fmpz_t res, ulong k, const fmpz_t n)

Sets `res` to $\sigma_k(|n|)$, the sum of the $k$-th powers of the
positive divisors of $n$, and to 0 for $n = 0$. Factors $n$.
:::

## Primality testing {#primality}

[](algorithms/number-theory/primality-proving) describes the proofs.

::: function
int fmpz_is_probabprime(const fmpz_t p)

Returns 0 if $p$ is composite and 1 if it is a probable prime. For $p$ up
to about 81 bits on 64-bit machines (where `n_is_prime` and
`n_ll_is_prime` decide primality), the answer is exact; otherwise the
test is trial division, a perfect square check and the
Baillie–PSW test (a base-2 strong probable prime test and a Lucas test),
for which no composite that passes is known.
:::

::: function
int fmpz_is_prime(const fmpz_t n)

Returns 1 if $n$ is prime and 0 otherwise, with a proof.

::: details
After the tests of `fmpz_is_probabprime` (which settle every $n$ up to
about 81 bits on 64-bit machines), the function trial factors $n - 1$ and $n + 1$ and tries the
Pocklington, Morrison and Brillhart–Lehmer–Selfridge tests [@BLS1975]. If
these do not succeed, the proof is done by ECPP (`ecpp_is_prime`) from 250
bits with one thread (the crossover rises with the number of threads)
and by APRCL (`aprcl_is_prime`) otherwise. If neither can decide, an exception is
raised; this is not expected to happen.
:::
:::

::: function
int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a)

Returns 1 if $n$ is a strong probable prime to base $a$, and 0 otherwise.
:::

::: function
int fmpz_is_probabprime_lucas(const fmpz_t n)

Returns 1 if $n$ is a Lucas probable prime with parameters chosen by
Selfridge's method A [@PSW1980], and 0 otherwise. No prime is declared
composite.
:::

::: function
int fmpz_is_probabprime_BPSW(const fmpz_t n)

Returns 1 if $n$ passes the Baillie–PSW test [@PSW1980] (a strong
probable prime test to base 2 and a Lucas test), and 0 otherwise. No
composite passing it is known, and no prime fails it.
:::

::: function
int fmpz_is_prime_pseudosquare(const fmpz_t n)

Returns 0 if $n$ is composite and 1 if it is proved prime by the
pseudosquare test (Theorem 2.7 of [@LukPatWil1996]), using a table of
pseudosquares up to $2^{64}$. Returns $-1$ if $n$ is too large (above
about 94 bits).
:::

::: function
int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, nn_ptr pm1, slong num_pm1)
int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, nn_ptr pp1, slong num_pp1)

The $n - 1$ test of Pocklington and the $n + 1$ test of Morrison for odd
$n$ (non-square for Morrison), given the `num_pm1` (`num_pp1`) primes
`pm1` (`pp1`) dividing $n - 1$ ($n + 1$). Return 0 if $n$ is composite.
Otherwise they set $F$ to the factored part of $n \mp 1$ and $R$ to
$(n \mp 1)/F$ and return 1, which proves that every prime factor of $n$ is
$\equiv 1 \pmod F$ (respectively $\equiv \pm 1$). This proves $n$ prime
if $F \ge \sqrt n$ (respectively $F > \sqrt n + 1$).
:::

::: function
void _fmpz_nm1_trial_factors(const fmpz_t n, nn_ptr pm1, slong * num_pm1, ulong limit)
void _fmpz_np1_trial_factors(const fmpz_t n, nn_ptr pp1, slong * num_pp1, ulong limit)

Trial factor $n - 1$ (respectively $n + 1$) by the primes up to about
`limit` and store the primes found in the array, setting its length. The
array needs room for $\log_2 n + 2$ entries. A `limit` of a small
multiple of $(\log n)^3 / 10$ is a reasonable choice for the tests above.
:::

::: function
void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t m, const fmpz_t n)
void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t m, const fmpz_t n)

Set $V_m$ and $V_{m+1}$ modulo $n$ for the Lucas sequence $V_0 = 2$,
$V_1 = A$, $V_j = A V_{j-1} - V_{j-2}$ (respectively $A V_{j-1} - B
V_{j-2}$), by the doubling formulas. No aliasing is allowed.
:::

::: function
void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t Un, const fmpz_t Un1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t k, const fmpz_t n)

For the Lucas sequence $U_0 = 0$, $U_1 = 1$, $U_j = A U_{j-1} - B
U_{j-2}$ modulo $n$: given $U_m, U_{m+1}$ (and $U_n, U_{n+1}$), compute
$U_{2m}, U_{2m+1}$, respectively $U_{m+n}, U_{m+n+1}$, respectively
$U_{km}, U_{km+1}$. Each output may alias the corresponding input ($U_{2m}$
with $U_m$, and so on); no other aliasing is allowed.
:::

::: function
void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, const fmpz_t Vm, const fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, const fmpz_t n)

Given $V_m$ and $V_{m+1}$ modulo $n$, sets $U_m$ and $U_{m+1}$, where
`Dinv` is the inverse of $D = A^2 - 4B$ modulo $n$. $U_m$ may alias $V_m$
and $U_{m+1}$ may alias $V_{m+1}$.
:::

::: function
void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved)

Sets `res` to the smallest prime greater than $n$ (2 for $n \le 1$). If
`proved` is nonzero, the result is proved prime with `fmpz_is_prime`;
otherwise it is a probable prime (`n_nextprime` for small $n$, GMP's
`mpz_nextprime`, a Baillie–PSW test, otherwise).
:::

::: function
int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s)

Lenstra's algorithm: for $0 < r < s < n$ with $\gcd(r, s) = 1$, finds a
divisor of $n$ congruent to $r$ modulo $s$ if there is one, sets `fac`
to it and returns 1; otherwise returns 0. Efficient when $s^3 > n$, as
in the final step of the $n \pm 1$ primality tests.
:::

## Internal functions {#internal}

These functions manage the representation directly. They are needed only
to write new low-level functions; see [](fmpz/representation). Those
returning `mpz_ptr` are declared only if `gmp.h` is included before the
FLINT headers.

::: function
mpz_ptr _fmpz_new_mpz(void)
void _fmpz_clear_mpz(fmpz f)

Take an initialised `mpz` with at least `MPZ_MIN_ALLOC` limbs from the
memory manager, respectively give back the `mpz` that the large value $f$
points to. The `mpz` given back must still have at least `MPZ_MIN_ALLOC`
limbs.
:::

::: function
void _fmpz_cleanup_mpz_content(void)
void _fmpz_cleanup(void)

Free the `mpz`s cached by the current thread (the second also frees the
cache array). Called by `flint_cleanup`; they do nothing in builds
configured with `--enable-reentrant`.
:::

::: function
mpz_ptr _fmpz_promote(fmpz_t f)
mpz_ptr _fmpz_promote_val(fmpz_t f)

Ensure that $f$ holds an `mpz` and return it. `_fmpz_promote` does not
preserve the value of a small $f$ (the `mpz` holds garbage); it is used
before overwriting $f$. `_fmpz_promote_val` preserves the value. Either
may leave $f$ temporarily non-canonical.
:::

::: function
void _fmpz_demote(fmpz_t f)
void _fmpz_demote_val(fmpz_t f)

`_fmpz_demote` releases the `mpz` of $f$, if any, and sets $f$ to 0 (any
value then fits inline). `_fmpz_demote_val` moves the value of $f$ back
inline if it fits, restoring canonical form.
:::

::: function
void _fmpz_promote_set_ui(fmpz_t f, ulong v)
void _fmpz_promote_set_si(fmpz_t f, slong v)
void _fmpz_promote_neg_ui(fmpz_t f, ulong v)
void _fmpz_promote_set_uiui(fmpz_t f, ulong hi, ulong lo)
void _fmpz_promote_neg_uiui(fmpz_t f, ulong hi, ulong lo)
void _fmpz_init_promote_set_ui(fmpz_t f, ulong v)
void _fmpz_init_promote_set_si(fmpz_t f, slong v)

Out-of-line parts of the inline assignment functions: set an initialised
(`_init_` versions: an uninitialised) $f$ to a value that does not fit
inline. Calling them with a value that fits leaves $f$ non-canonical.
:::

::: function
void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z)
void _fmpz_clear_readonly_mpz(mpz_t z)

Helpers of the read-only conversions: wrap the `mpz` $z$ in a fresh `mpz`
slot of $f$, respectively release the `mpz` of such an $f$ when its value
is small.
:::

::: function
int _fmpz_is_canonical(const fmpz_t x)

Returns whether $x$ is in canonical form: inline, or an `mpz` whose
value does not fit inline and whose top limb is nonzero. Used in tests.
:::

::: function
void _fmpz_smod(fmpz_t r, const fmpz_t a, const fmpz_t m, int sign, fmpz_t t)

Sets $r$ to $a$ modulo $m \ne 0$ in a range selected by `sign`: for
`sign` $< 0$ the remainder of truncating division ($|r| < |m|$, sign of
$a$), for 0 that of floor division ($|r| < |m|$, sign of $m$), and for
`sign` $> 0$ the symmetric remainder $-|m| < 2r \le |m|$. The scratch
variable $t$ must not alias any other argument.
:::

::: function
slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv)

Divides $x \ne 0$ in place by the largest power of $f > 1$ dividing it
and returns the exponent. `finv` must be `n_precompute_inverse(f)` if $f$
is small and 0 otherwise. $x$ and $f$ must be distinct.
:::
