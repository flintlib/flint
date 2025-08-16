.. _arf:

**arf.h** -- arbitrary-precision floating-point numbers
===============================================================================

A variable of type :type:`arf_t` holds an arbitrary-precision binary
floating-point number: that is, a rational number of the form
`x \cdot 2^y` where `x, y \in \mathbb{Z}` and `x` is odd,
or one of the special values zero, plus infinity, minus infinity,
or NaN (not-a-number).
There is currently no support for negative zero, unsigned infinity,
or a NaN with a payload.

The *exponent* of a finite and nonzero floating-point number can be
defined in different
ways: for example, as the component *y* above, or as the unique
integer *e* such that
`x \cdot 2^y = m \cdot 2^e` where `0.5 \le |m| < 1`.
The internal representation of an :type:`arf_t` stores the
exponent in the latter format.

Except where otherwise noted, functions have the following semantics:

* Functions taking *prec* and *rnd* parameters at the end of the argument list
  and returning an ``int`` flag round the result in the output variable
  to *prec* bits in the direction specified by *rnd*. The return flag
  is 0 if the result is exact
  (not rounded) and 1 if the result is inexact (rounded).
  Correct rounding is guaranteed: the result is the floating-point number
  obtained by viewing the inputs as exact numbers, in principle carrying out
  the mathematical operation exactly, and rounding the resulting real number
  to the nearest representable floating-point number whose mantissa has at
  most the specified number of bits, in the specified direction of rounding.
  In particular, the error is at most 1 ulp with directed rounding modes
  and 0.5 ulp when rounding to nearest.

* Other functions perform the operation exactly.

Since exponents are bignums, overflow or underflow cannot occur.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arf_struct

.. type:: arf_t

    An :type:`arf_struct` contains four words: an :type:`fmpz` exponent (*exp*),
    a *size* field tracking the number of limbs used (one bit of this
    field is also used for the sign of the number), and two more words.
    The last two words hold the value directly if there are at most two limbs,
    and otherwise contain one *alloc* field (tracking the total number of
    allocated limbs, not all of which might be used) and a pointer to
    the actual limbs.
    Thus, up to 128 bits on a 64-bit machine and 64 bits on a 32-bit machine,
    no space outside of the :type:`arf_struct` is used.

    An :type:`arf_t` is defined as an array of length one of type
    :type:`arf_struct`, permitting an :type:`arf_t` to be passed by reference.

.. type:: arf_rnd_t

    Specifies the rounding mode for the result of an approximate operation.

.. macro:: ARF_RND_DOWN

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards zero.

.. macro:: ARF_RND_UP

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction away from zero.

.. macro:: ARF_RND_FLOOR

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards minus infinity.

.. macro:: ARF_RND_CEIL

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards plus infinity.

.. macro:: ARF_RND_NEAR

    Specifies that the result of an operation should be rounded to the
    nearest representable number, rounding to even if there is a tie
    between two values.

.. macro:: ARF_PREC_EXACT

    If passed as the precision parameter to a function, indicates that no
    rounding is to be performed. **Warning**: use of this value is unsafe in
    general. It must only be
    passed as input under the following two conditions:

    * The operation in question can inherently be viewed as an exact operation
      in `\mathbb{Z}[\tfrac{1}{2}]` for all possible inputs, provided that
      the precision is large enough. Examples include addition,
      multiplication, conversion from integer types to arbitrary-precision
      floating-point types, and evaluation of some integer-valued functions.

    * The exact result of the operation will certainly fit in memory.
      Note that, for example, adding two numbers whose exponents are far
      apart can easily produce an exact result that is far too large to
      store in memory.

    The typical use case is to work with small integer values, double
    precision constants, and the like. It is also useful when writing
    test code. If in doubt, simply try with some convenient high precision
    instead of using this special value, and check that the result is exact.

Memory management
-------------------------------------------------------------------------------

.. function:: void arf_init(arf_t x)

    Initializes the variable *x* for use. Its value is set to zero.

.. function:: void arf_clear(arf_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: slong arf_allocated_bytes(const arf_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(arf_struct)`` to get the size of the object as a whole.

Special values
-------------------------------------------------------------------------------

.. function:: void arf_zero(arf_t res)

.. function:: void arf_one(arf_t res)

.. function:: void arf_pos_inf(arf_t res)

.. function:: void arf_neg_inf(arf_t res)

.. function:: void arf_nan(arf_t res)

    Sets *res* respectively to 0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_zero(const arf_t x)

.. function:: int arf_is_one(const arf_t x)

.. function:: int arf_is_pos_inf(const arf_t x)

.. function:: int arf_is_neg_inf(const arf_t x)

.. function:: int arf_is_nan(const arf_t x)

    Returns nonzero iff *x* respectively equals 0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_inf(const arf_t x)

    Returns nonzero iff *x* equals either `+\infty` or `-\infty`.

.. function:: int arf_is_normal(const arf_t x)

    Returns nonzero iff *x* is a finite, nonzero floating-point value, i.e.
    not one of the special values 0, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_special(const arf_t x)

    Returns nonzero iff *x* is one of the special values
    0, `+\infty`, `-\infty`, NaN, i.e. not a finite, nonzero
    floating-point value.

.. function:: int arf_is_finite(const arf_t x)

    Returns nonzero iff *x* is a finite floating-point value,
    i.e. not one of the values `+\infty`, `-\infty`, NaN.
    (Note that this is not equivalent to the negation of
    :func:`arf_is_inf`.)


Assignment, rounding and conversions
-------------------------------------------------------------------------------

.. function:: void arf_set(arf_t res, const arf_t x)

.. function:: void arf_set_mpz(arf_t res, const mpz_t x)

.. function:: void arf_set_fmpz(arf_t res, const fmpz_t x)

.. function:: void arf_set_ui(arf_t res, ulong x)

.. function:: void arf_set_si(arf_t res, slong x)

.. function:: void arf_set_mpfr(arf_t res, const mpfr_t x)

.. function:: void arf_set_d(arf_t res, double x)

    Sets *res* to the exact value of *x*.

.. function:: void arf_swap(arf_t x, arf_t y)

    Swaps *x* and *y* efficiently.

.. function:: void arf_init_set_ui(arf_t res, ulong x)

.. function:: void arf_init_set_si(arf_t res, slong x)

    Initializes *res* and sets it to *x* in a single operation.

.. function:: int arf_set_round(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

.. function:: int arf_set_round_si(arf_t res, slong x, slong prec, arf_rnd_t rnd)

.. function:: int arf_set_round_ui(arf_t res, ulong x, slong prec, arf_rnd_t rnd)

.. function:: int arf_set_round_mpz(arf_t res, const mpz_t x, slong prec, arf_rnd_t rnd)

.. function:: int arf_set_round_fmpz(arf_t res, const fmpz_t x, slong prec, arf_rnd_t rnd)

    Sets *res* to *x*, rounded to *prec* bits in the direction
    specified by *rnd*.

.. function:: void arf_set_si_2exp_si(arf_t res, slong m, slong e)

.. function:: void arf_set_ui_2exp_si(arf_t res, ulong m, slong e)

.. function:: void arf_set_fmpz_2exp(arf_t res, const fmpz_t m, const fmpz_t e)

    Sets *res* to `m \cdot 2^e`.

.. function:: int arf_set_round_fmpz_2exp(arf_t res, const fmpz_t x, const fmpz_t e, slong prec, arf_rnd_t rnd)

    Sets *res* to `x \cdot 2^e`, rounded to *prec* bits in the direction
    specified by *rnd*.

.. function:: void arf_get_fmpz_2exp(fmpz_t m, fmpz_t e, const arf_t x)

    Sets *m* and *e* to the unique integers such that
    `x = m \cdot 2^e` and *m* is odd,
    provided that *x* is a nonzero finite fraction.
    If *x* is zero, both *m* and *e* are set to zero. If *x* is
    infinite or NaN, the result is undefined.

.. function:: void arf_frexp(arf_t m, fmpz_t e, const arf_t x)

    Writes *x* as `m \cdot 2^e`, where `0.5 \le |m| < 1` if *x* is a normal
    value. If *x* is a special value, copies this to *m* and sets *e* to zero.
    Note: for the inverse operation (*ldexp*), use :func:`arf_mul_2exp_fmpz`.

.. function:: double arf_get_d(const arf_t x, arf_rnd_t rnd)

    Returns *x* rounded to a double in the direction specified by *rnd*.
    This method rounds correctly when overflowing or underflowing
    the double exponent range (this was not the case in an earlier version).

.. function:: int arf_get_mpfr(mpfr_t res, const arf_t x, mpfr_rnd_t rnd)

    Sets the MPFR variable *res* to the value of *x*. If the precision of *x*
    is too small to allow *res* to be represented exactly, it is rounded in
    the specified MPFR rounding mode. The return value (-1, 0 or 1)
    indicates the direction of rounding, following the convention
    of the MPFR library.

    If *x* has an exponent too large or small to fit in the MPFR type, the
    result overflows to an infinity or underflows to a (signed) zero,
    and the corresponding MPFR exception flags are set.

.. function:: int arf_get_fmpz(fmpz_t res, const arf_t x, arf_rnd_t rnd)

    Sets *res* to *x* rounded to the nearest integer in the direction
    specified by *rnd*. If rnd is *ARF_RND_NEAR*, rounds to the nearest
    even integer in case of a tie. Returns inexact (beware: accordingly
    returns whether *x* is *not* an integer).

    This method aborts if *x* is infinite or NaN, or if the exponent of *x*
    is so large that allocating memory for the result fails.

    Warning: this method will allocate a huge amount of memory to store
    the result if the exponent of *x* is huge. Memory allocation could
    succeed even if the required space is far larger than the physical
    memory available on the machine, resulting in swapping. It is recommended
    to check that *x* is within a reasonable range before calling this method.

.. function:: slong arf_get_si(const arf_t x, arf_rnd_t rnd)

    Returns *x* rounded to the nearest integer in the direction specified by
    *rnd*. If *rnd* is *ARF_RND_NEAR*, rounds to the nearest even integer
    in case of a tie. Aborts if *x* is infinite, NaN, or the value is
    too large to fit in a slong.

.. function:: int arf_get_fmpz_fixed_fmpz(fmpz_t res, const arf_t x, const fmpz_t e)

.. function:: int arf_get_fmpz_fixed_si(fmpz_t res, const arf_t x, slong e)

    Converts *x* to a mantissa with predetermined exponent, i.e. sets *res* to
    an integer *y* such that `y \times 2^e \approx x`, truncating if necessary.
    Returns 0 if exact and 1 if truncation occurred.

    The warnings for :func:`arf_get_fmpz` apply.

.. function:: void arf_floor(arf_t res, const arf_t x)
              void arf_ceil(arf_t res, const arf_t x)
              void arf_nint(arf_t res, const arf_t x)

    Sets *res* respectively to the floor `\lfloor x \rfloor`,
    ceiling `\lceil x \rceil`
    or nearest integer `\lfloor x \rceil` (with ties-to-even).
    The result is always represented exactly, requiring no more bits to
    store than the input. To round the result to a floating-point number
    with a lower precision, call :func:`arf_set_round` afterwards.

.. function:: void arf_get_fmpq(fmpq_t res, const arf_t x)

    Set *res* to the exact rational value of *x*.
    This method aborts if *x* is infinite or NaN, or if the exponent of *x*
    is so large that allocating memory for the result fails.

Comparisons and bounds
-------------------------------------------------------------------------------

.. function:: int arf_equal(const arf_t x, const arf_t y)
              int arf_equal_si(const arf_t x, slong y)
              int arf_equal_ui(const arf_t x, ulong y)
              int arf_equal_d(const arf_t x, double y)

    Returns nonzero iff *x* and *y* are exactly equal. NaN is not
    treated specially, i.e. NaN compares as equal to itself.

    For comparison with a *double*, the values -0 and +0 are
    both treated as zero, and all NaN values are treated as identical.

.. function:: int arf_cmp(const arf_t x, const arf_t y)

.. function:: int arf_cmp_si(const arf_t x, slong y)

.. function:: int arf_cmp_ui(const arf_t x, ulong y)

.. function:: int arf_cmp_d(const arf_t x, double y)

    Returns negative, zero, or positive, depending on whether *x* is
    respectively smaller, equal, or greater compared to *y*.
    Comparison with NaN is undefined.

.. function:: int arf_cmpabs(const arf_t x, const arf_t y)

.. function:: int arf_cmpabs_ui(const arf_t x, ulong y)

.. function:: int arf_cmpabs_d(const arf_t x, double y)

.. function:: int arf_cmpabs_mag(const arf_t x, const mag_t y)

    Compares the absolute values of *x* and *y*.

.. function:: int arf_cmp_2exp_si(const arf_t x, slong e)

.. function:: int arf_cmpabs_2exp_si(const arf_t x, slong e)

    Compares *x* (respectively its absolute value) with `2^e`.

.. function:: int arf_sgn(const arf_t x)

    Returns `-1`, `0` or `+1` according to the sign of *x*. The sign
    of NaN is undefined.

.. function:: void arf_min(arf_t res, const arf_t a, const arf_t b)

.. function:: void arf_max(arf_t res, const arf_t a, const arf_t b)

    Sets *res* respectively to the minimum and the maximum of *a* and *b*.

.. function:: slong arf_bits(const arf_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of *x*, i.e. the minimum precision sufficient to represent
    *x* exactly. Returns 0 if *x* is a special value.

.. function:: int arf_is_int(const arf_t x)

    Returns nonzero iff *x* is integer-valued.

.. function:: int arf_is_int_2exp_si(const arf_t x, slong e)

    Returns nonzero iff *x* equals `n 2^e` for some integer *n*.

.. function:: void arf_abs_bound_lt_2exp_fmpz(fmpz_t res, const arf_t x)

    Sets *res* to the smallest integer *b* such that `|x| < 2^b`.
    If *x* is zero, infinity or NaN, the result is undefined.

.. function:: void arf_abs_bound_le_2exp_fmpz(fmpz_t res, const arf_t x)

    Sets *res* to the smallest integer *b* such that `|x| \le 2^b`.
    If *x* is zero, infinity or NaN, the result is undefined.

.. function:: slong arf_abs_bound_lt_2exp_si(const arf_t x)

    Returns the smallest integer *b* such that `|x| < 2^b`, clamping
    the result to lie between -*ARF_PREC_EXACT* and *ARF_PREC_EXACT*
    inclusive. If *x* is zero, -*ARF_PREC_EXACT* is returned,
    and if *x* is infinity or NaN, *ARF_PREC_EXACT* is returned.

Magnitude functions
-------------------------------------------------------------------------------

.. function:: void arf_get_mag(mag_t res, const arf_t x)

    Sets *res* to an upper bound for the absolute value of *x*.

.. function:: void arf_get_mag_lower(mag_t res, const arf_t x)

    Sets *res* to a lower bound for the absolute value of *x*.

.. function:: void arf_set_mag(arf_t res, const mag_t x)

    Sets *res* to *x*. This operation is exact.

.. function:: void mag_init_set_arf(mag_t res, const arf_t x)

    Initializes *res* and sets it to an upper bound for *x*.

.. function:: void mag_fast_init_set_arf(mag_t res, const arf_t x)

    Initializes *res* and sets it to an upper bound for *x*.
    Assumes that the exponent of *res* is small (this function is unsafe).

.. function:: void arf_mag_set_ulp(mag_t res, const arf_t x, slong prec)

    Sets *res* to the magnitude of the unit in the last place (ulp) of *x*
    at precision *prec*.

.. function:: void arf_mag_add_ulp(mag_t res, const mag_t x, const arf_t y, slong prec)

    Sets *res* to an upper bound for the sum of *x* and the
    magnitude of the unit in the last place (ulp) of *y*
    at precision *prec*.

.. function:: void arf_mag_fast_add_ulp(mag_t res, const mag_t x, const arf_t y, slong prec)

    Sets *res* to an upper bound for the sum of *x* and the
    magnitude of the unit in the last place (ulp) of *y*
    at precision *prec*. Assumes that all exponents are small.

Shallow assignment
-------------------------------------------------------------------------------

.. function:: void arf_init_set_shallow(arf_t z, const arf_t x)

.. function:: void arf_init_set_mag_shallow(arf_t z, const mag_t x)

    Initializes *z* to a shallow copy of *x*. A shallow copy just involves
    copying struct data (no heap allocation is performed).

    The target variable *z* may not be cleared or modified in any way (it can
    only be used as constant input to functions), and may not be used after
    *x* has been cleared. Moreover, after *x* has been assigned shallowly
    to *z*, no modification of *x* is permitted as slong as *z* is in use.

.. function:: void arf_init_neg_shallow(arf_t z, const arf_t x)

.. function:: void arf_init_neg_mag_shallow(arf_t z, const mag_t x)

    Initializes *z* shallowly to the negation of *x*.

Random number generation
-------------------------------------------------------------------------------

.. function:: void arf_randtest(arf_t res, flint_rand_t state, slong bits, slong mag_bits)

    Generates a finite random number whose mantissa has precision at most
    *bits* and whose exponent has at most *mag_bits* bits. The
    values are distributed non-uniformly: special bit patterns are generated
    with high probability in order to allow the test code to exercise corner
    cases.

.. function:: void arf_randtest_not_zero(arf_t res, flint_rand_t state, slong bits, slong mag_bits)

    Identical to :func:`arf_randtest`, except that zero is never produced
    as an output.

.. function:: void arf_randtest_special(arf_t res, flint_rand_t state, slong bits, slong mag_bits)

    Identical to :func:`arf_randtest`, except that the output occasionally
    is set to an infinity or NaN.

.. function:: void arf_urandom(arf_t res, flint_rand_t state, slong bits, arf_rnd_t rnd)

    Sets *res* to a uniformly distributed random number in the interval
    `[0, 1]`. The method uses rounding from integers to floats based on the
    rounding mode *rnd*.

Input and output
-------------------------------------------------------------------------------

.. function:: void arf_debug(const arf_t x)

    Prints information about the internal representation of *x*.

.. function:: void arf_print(const arf_t x)

    Prints *x* as an integer mantissa and exponent.

.. function:: void arf_printd(const arf_t x, slong d)

    Prints *x* as a decimal floating-point number, rounding to *d* digits.
    Rounding is faithful (at most 1 ulp error).

.. function:: char * arf_get_str(const arf_t x, slong d)

    Returns *x* as a decimal floating-point number, rounding to *d* digits.
    Rounding is faithful (at most 1 ulp error).

.. function:: void arf_fprint(FILE * file, const arf_t x)

    Prints *x* as an integer mantissa and exponent to the stream *file*.

.. function:: void arf_fprintd(FILE * file, const arf_t y, slong d)

    Prints *x* as a decimal floating-point number to the stream *file*,
    rounding to *d* digits.
    Rounding is faithful (at most 1 ulp error).

.. function:: char * arf_dump_str(const arf_t x)

    Allocates a string and writes a binary representation of *x* to it that can
    be read by :func:`arf_load_str`. The returned string needs to be
    deallocated with *flint_free*.

.. function:: int arf_load_str(arf_t x, const char * str)

    Parses *str* into *x*. Returns a nonzero value if *str* is not formatted
    correctly.

.. function:: int arf_dump_file(FILE * stream, const arf_t x)

    Writes a binary representation of *x* to *stream* that can be read by
    :func:`arf_load_file`. Returns a nonzero value if the data could not be
    written.

.. function:: int arf_load_file(arf_t x, FILE * stream)

    Reads *x* from *stream*. Returns a nonzero value if the data is not
    formatted correctly or the read failed. Note that the data is assumed to be
    delimited by a whitespace or end-of-file, i.e., when writing multiple
    values with :func:`arf_dump_file` make sure to insert a whitespace to
    separate consecutive values.

Addition and multiplication
-------------------------------------------------------------------------------

.. function:: void arf_abs(arf_t res, const arf_t x)

    Sets *res* to the absolute value of *x* exactly.

.. function:: void arf_neg(arf_t res, const arf_t x)

    Sets *res* to `-x` exactly.

.. function:: int arf_neg_round(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

    Sets *res* to `-x`.

.. function:: int arf_add(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_add_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_add_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_add_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Sets *res* to `x + y`.

.. function:: int arf_add_fmpz_2exp(arf_t res, const arf_t x, const fmpz_t y, const fmpz_t e, slong prec, arf_rnd_t rnd)

    Sets *res* to `x + y 2^e`.

.. function:: int arf_sub(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_sub_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_sub_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_sub_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Sets *res* to `x - y`.

.. function:: void arf_mul_2exp_si(arf_t res, const arf_t x, slong e)

.. function:: void arf_mul_2exp_fmpz(arf_t res, const arf_t x, const fmpz_t e)

    Sets *res* to `x 2^e` exactly.

.. function:: int arf_mul(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_mul_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_mul_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_mul_mpz(arf_t res, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_mul_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Sets *res* to `x \cdot y`.

.. function:: int arf_addmul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_addmul_ui(arf_t z, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_addmul_si(arf_t z, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_addmul_mpz(arf_t z, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_addmul_fmpz(arf_t z, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Performs a fused multiply-add `z = z + x \cdot y`, updating *z* in-place.

.. function:: int arf_submul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_submul_ui(arf_t z, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_submul_si(arf_t z, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_submul_mpz(arf_t z, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_submul_fmpz(arf_t z, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Performs a fused multiply-subtract `z = z - x \cdot y`, updating *z* in-place.

.. function:: int arf_fma(arf_t res, const arf_t x, const arf_t y, const arf_t z, slong prec, arf_rnd_t rnd)

    Sets *res* to `x \cdot y + z`. This is equivalent to an *addmul* except
    that *res* and *z* can be separate variables.

.. function:: int arf_sosq(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    Sets *res* to `x^2 + y^2`, rounded to *prec* bits in the direction specified by *rnd*.

Summation
-------------------------------------------------------------------------------

.. function:: int arf_sum(arf_t res, arf_srcptr terms, slong len, slong prec, arf_rnd_t rnd)

    Sets *res* to the sum of the array *terms* of length *len*, rounded to
    *prec* bits in the direction specified by *rnd*. The sum is computed as if
    done without any intermediate rounding error, with only a single rounding
    applied to the final result. Unlike repeated calls to :func:`arf_add` with
    infinite precision, this function does not overflow if the magnitudes of
    the terms are far apart. Warning: this function is implemented naively,
    and the running time is quadratic with respect to *len* in the worst case.

Dot products
-------------------------------------------------------------------------------

.. function:: int arf_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

    Computes a dot product, with the same meaning of
    the parameters as :func:`arb_dot`. The result is correctly rounded
    in the direction specified by ``rnd``. This function is currently not optimized
    for speed.

.. function:: void arf_approx_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

    Computes an approximate dot product, with the same meaning of
    the parameters as :func:`arb_dot`.
    This operation is not correctly rounded: the final rounding is done
    in the direction ``rnd`` but intermediate roundings are
    implementation-defined.

Division
-------------------------------------------------------------------------------

.. function:: int arf_div(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_div_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_ui_div(arf_t res, ulong x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_div_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

.. function:: int arf_si_div(arf_t res, slong x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_div_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_fmpz_div(arf_t res, const fmpz_t x, const arf_t y, slong prec, arf_rnd_t rnd)

.. function:: int arf_fmpz_div_fmpz(arf_t res, const fmpz_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    Sets *res* to `x / y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact. The result is NaN if *y* is zero.

Square roots
-------------------------------------------------------------------------------

.. function:: int arf_sqrt(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

.. function:: int arf_sqrt_ui(arf_t res, ulong x, slong prec, arf_rnd_t rnd)

.. function:: int arf_sqrt_fmpz(arf_t res, const fmpz_t x, slong prec, arf_rnd_t rnd)

    Sets *res* to `\sqrt{x}`. The result is NaN if *x* is negative.

.. function:: int arf_rsqrt(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

    Sets *res* to `1/\sqrt{x}`. The result is NaN if *x* is
    negative, and `+\infty` if *x* is zero.

.. function:: int arf_root(arf_t res, const arf_t x, ulong k, slong prec, arf_rnd_t rnd)

    Sets *res* to `x^{1/k}`. The result is NaN if *x* is negative.
    Warning: this function is a wrapper around the MPFR root function.
    It gets slow and uses much memory for large *k*.
    Consider working with :func:`arb_root_ui` for large *k* instead of using this
    function directly.

Complex arithmetic
-------------------------------------------------------------------------------

.. function:: int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)

.. function:: int arf_complex_mul_fallback(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)

    Computes the complex product `e + fi = (a + bi)(c + di)`, rounding both
    `e` and `f` correctly to *prec* bits in the direction specified by *rnd*.
    The first bit in the return code indicates inexactness of `e`, and the
    second bit indicates inexactness of `f`.

    If any of the components *a*, *b*, *c*, *d* is zero, two real
    multiplications and no additions are done. This convention is used even
    if any other part contains an infinity or NaN, and the behavior
    with infinite/NaN input is defined accordingly.

    The *fallback* version is implemented naively, for testing purposes.
    No squaring optimization is implemented.

.. function:: int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)

    Computes the complex square `e + fi = (a + bi)^2`. This function has
    identical semantics to :func:`arf_complex_mul` (with `c = a, b = d`),
    but is faster.

Low-level methods
-------------------------------------------------------------------------------

.. function:: int _arf_get_integer_mpn(nn_ptr y, nn_srcptr xp, slong xn, slong exp)

    Given a floating-point number *x* represented by *xn* limbs at *xp*
    and an exponent *exp*, writes the integer part of *x* to
    *y*, returning whether the result is inexact.
    The correct number of limbs is written (no limbs are written
    if the integer part of *x* is zero).
    Assumes that ``xp[0]`` is nonzero and that the
    top bit of ``xp[xn-1]`` is set.

.. function:: int _arf_set_mpn_fixed(arf_t z, nn_srcptr xp, slong xn, slong fixn, int negative, slong prec, arf_rnd_t rnd)

    Sets *z* to the fixed-point number having *xn* total limbs and *fixn*
    fractional limbs, negated if *negative* is set, rounding *z* to *prec*
    bits in the direction *rnd* and returning whether the result is inexact.
    Both *xn* and *fixn* must be nonnegative and not so large
    that the bit shift would overflow an *slong*, but otherwise no
    assumptions are made about the input.

.. function:: int _arf_set_round_ui(arf_t z, ulong x, int sgnbit, slong prec, arf_rnd_t rnd)

    Sets *z* to the integer *x*, negated if *sgnbit* is 1, rounded to *prec*
    bits in the direction specified by *rnd*. There are no assumptions on *x*.

.. function:: int _arf_set_round_uiui(arf_t z, slong * fix, ulong hi, ulong lo, int sgnbit, slong prec, arf_rnd_t rnd)

    Sets the mantissa of *z* to the two-limb mantissa given by *hi* and *lo*,
    negated if *sgnbit* is 1, rounded to *prec* bits in the direction specified
    by *rnd*. Requires that not both *hi* and *lo* are zero.
    Writes the exponent shift to *fix* without writing the exponent of *z*
    directly.

.. function:: int _arf_set_round_mpn(arf_t z, slong * exp_shift, nn_srcptr x, slong xn, int sgnbit, slong prec, arf_rnd_t rnd)

    Sets the mantissa of *z* to the mantissa given by the *xn* limbs in *x*,
    negated if *sgnbit* is 1, rounded to *prec* bits in the direction
    specified by *rnd*. Returns the inexact flag. Requires that *xn* is positive
    and that the top limb of *x* is nonzero. If *x* has leading zero bits,
    writes the shift to *exp_shift*. This method does not write the exponent of
    *z* directly. Requires that *x* does not point to the limbs of *z*.

