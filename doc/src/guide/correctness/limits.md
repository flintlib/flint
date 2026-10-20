---
title: Sizes and limits
summary: Why lengths, shift counts, exponents and precisions are not checked for overflow, examples of what goes wrong with huge or negative values, the size limits inherited from GMP, exponent ranges, and what happens when memory runs out.
order: 4
---

FLINT does not check size parameters (lengths of arrays and polynomials,
dimensions, shift counts, exponents, precisions) for values that are
negative, absurdly large, or large enough to overflow the arithmetic
that computes memory sizes. Such values give undefined behaviour: a
failed allocation, a crash, or a silently wrong result. On a 64-bit
machine a legitimate computation runs out of memory or time long before
it reaches these limits, so the problem arises mainly with sizes that
come from a programming error (an uninitialised variable, a negative
value passed as a length, a runaway loop that keeps doubling a
precision) or from untrusted input.

## Overflow in size computations {#overflow}

Lengths, degrees, dimensions and precisions are `slong` values, shift
counts and exponents usually `ulong` or `flint_bitcnt_t`, and the
number of bytes to allocate is computed from them in word arithmetic,
which wraps around silently. Two examples on a 64-bit machine:

- `nmod_poly_init2(f, n, WORD(1) << 61)` asks for $2^{61}$
  coefficients of 8 bytes. The byte count $2^{64}$ wraps around to 0,
  the allocation succeeds, and the polynomial claims an allocation of
  $2^{61}$ coefficients that does not exist: the first write to it
  corrupts the heap.
- `fmpz_mul_2exp(x, y, e)` computes the bit length of the result as
  `fmpz_bits(y) + e` in a `ulong`. For $y = 1$ and $e = 2^{64} - 16$ this
  wraps around, and the function silently returns a 1-bit number. For
  $e$ slightly above $2^{37}$ the number of limbs, $e/64$, no longer
  fits in the `int` that GMP uses for sizes, and the function crashes.
  For $e = 2^{36}$ it tries to allocate 8 GB, which on a machine with
  less memory ends in an allocation failure.

The same applies to negative values: a negative length passed to a
function that expects a nonnegative one is converted to a huge
unsigned size, or used as a loop bound that is never reached.

## Size limits of integers {#integers}

FLINT's integers (`fmpz_t`) are GMP integers when they do not fit in a
word, and GMP stores the number of limbs of an integer in an `int`. An
integer therefore has at most $2^{31} - 1$ limbs: about $2^{37}$ bits,
or 41 billion decimal digits, on a 64-bit machine. Some FLINT functions
compute limb counts in `int` themselves and fail somewhat earlier, as
`fmpz_mul_2exp` above. Everything built on `fmpz_t` inherits the limit:
the coefficients of `fmpz_poly`, the numerators of `fmpq`, the mantissas
of large `arf_t` numbers, the precision of computations that convert to
GMP integers internally. In practice the limit only matters on machines
with hundreds of gigabytes of memory, but it is not checked.

On 32-bit machines the address space, not GMP, is the limit, and sizes
are much smaller: precisions and lengths above $2^{24}$ or so should be
regarded as unsupported, since intermediate size computations in
`slong` can overflow.

## Exponents {#exponents}

Several types have exponents of unbounded size, so that overflow cannot
occur:

- `arf_t`, `mag_t` and therefore `arb_t` and `acb_t` store their
  exponents as `fmpz` integers, so a ball such as $[2^{2^{100}} \pm 1]$
  is cheap to represent. Conversions whose output size depends on the
  exponent, such as `arf_get_fmpz`, are not protected; see
  [](guide/types/ball-arithmetic#limits).
- The multivariate polynomial types (`fmpz_mpoly` and the others) pack
  exponents into fields of a variable number of bits and repack
  automatically when an operation needs larger fields, up to
  multiprecision exponents. Some conversions to types with word-size
  exponents check and throw, for example `fmpz_mpoly_to_fmpz_poly`
  raises `FLINT_EXPOF` when the exponents of its input use more than
  one word.

Other types have bounded exponents and report overflow: the `nfloat`
types have exponents in $\pm$`WORD_MAX / 4`, and an operation whose
result overflows returns `GR_UNABLE` (or an infinity, if the context
allows infinities). Word-size types such as `ulong` and `double` have
their usual limits, and the functions that take or return them
(`fmpz_get_si`, `arf_get_d`, ...) document what happens out of range.

## Running out of memory {#memory}

All memory that FLINT allocates itself comes from `flint_malloc`,
`flint_calloc`, `flint_realloc` and `flint_aligned_alloc`. These never
return `NULL`: when the allocation fails they call `flint_throw` with the
message "Unable to allocate memory", which by default prints it and
aborts the program. The memory of large `fmpz` integers is allocated by
GMP, whose default allocation functions print "GNU MP: Cannot allocate
memory" (or "reallocate") to standard error and abort. A program that
installs an exception handler with `flint_set_throw` and jumps out of
it with `longjmp` should not continue to use the objects involved in
the failed operation, and memory allocated by that operation is lost.

On systems that overcommit memory, such as Linux in its default
configuration, a huge allocation may "succeed" and the process be
killed later, when the memory is first used. To make sure that
allocations fail early and visibly, limit the address space of the
process (with `ulimit -v` in the shell or `setrlimit` in C).

Stack space is limited too. Functions take small temporary blocks
from the stack (`TMP_ALLOC` takes blocks of up to 8192 bytes there);
with the usual stack sizes this is not a problem, but threads created
with a very small stack can overflow it.

## Practical advice {#advice}

- Compute sizes in a type that cannot overflow, or check them before
  passing them to FLINT: a length should be nonnegative and at most
  `WORD_MAX` divided by the size of an element.
- In loops that increase a precision or a length, put an upper limit
  on the value instead of looping until success.
- When sizes come from outside the program, check them against
  explicit bounds that fit the available memory and time; see
  [](guide/correctness/security#untrusted).
