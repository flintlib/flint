.. _fmpq:

**fmpq.h** -- rational numbers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpq

.. type:: fmpq_t

    Description.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpq_init(fmpq_t x)

    Initialises the ``fmpq_t`` variable ``x`` for use. Its value
    is set to 0.

.. function:: void fmpq_clear(fmpq_t x)

    Clears the ``fmpq_t`` variable ``x``. To use the variable again,
    it must be re-initialised with ``fmpq_init``.


Canonicalisation
--------------------------------------------------------------------------------


.. function:: void fmpq_canonicalise(fmpq_t res)

    Puts ``res`` in canonical form: the numerator and denominator are
    reduced to lowest terms, and the denominator is made positive.
    If the numerator is zero, the denominator is set to one.

    If the denominator is zero, the outcome of calling this function is
    undefined, regardless of the value of the numerator.

.. function:: void _fmpq_canonicalise(fmpz_t num, fmpz_t den)

    Does the same thing as ``fmpq_canonicalise``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    of ``num`` and ``den`` is not allowed.

.. function:: int fmpq_is_canonical(const fmpq_t x)

    Returns nonzero if ``fmpq_t`` x is in canonical form
    (as produced by ``fmpq_canonicalise``), and zero otherwise.

.. function:: int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)

    Does the same thing as ``fmpq_is_canonical``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables.


Basic assignment
--------------------------------------------------------------------------------


.. function:: void fmpq_set(fmpq_t dest, const fmpq_t src)

    Sets ``dest`` to a copy of ``src``. No canonicalisation
    is performed.

.. function:: void fmpq_swap(fmpq_t op1, fmpq_t op2)

    Swaps the two rational numbers ``op1`` and ``op2``.

.. function:: void fmpq_neg(fmpq_t dest, const fmpq_t src)

    Sets ``dest`` to the additive inverse of ``src``.

.. function:: void fmpq_abs(fmpq_t dest, const fmpq_t src)

    Sets ``dest`` to the absolute value of ``src``.

.. function:: void fmpq_zero(fmpq_t res)

    Sets the value of ``res`` to 0.

.. function:: void fmpq_one(fmpq_t res)

    Sets the value of ``res`` to `1`.


Comparison
--------------------------------------------------------------------------------


.. function:: int fmpq_is_zero(const fmpq_t res)

    Returns nonzero if ``res`` has value 0, and returns zero otherwise.

.. function:: int fmpq_is_one(const fmpq_t res)

    Returns nonzero if ``res`` has value `1`, and returns zero otherwise.

.. function:: int fmpq_is_pm1(const fmpq_t res)

    Returns nonzero if ``res`` has value `\pm{1}` and zero otherwise.

.. function:: int fmpq_equal(const fmpq_t x, const fmpq_t y)

    Returns nonzero if ``x`` and ``y`` are equal, and zero otherwise.
    Assumes that ``x`` and ``y`` are both in canonical form.

.. function:: int fmpq_sgn(const fmpq_t x)

    Returns the sign of the rational number `x`.

.. function:: int fmpq_cmp(const fmpq_t x, const fmpq_t y)

    Returns negative if `x < y`, zero if `x = y`, and positive if `x > y`.

.. function:: int fmpq_cmp_ui(const fmpq_t x, ulong y)

    Returns negative if `x < y`, zero if `x = y`, and positive if `x > y`.

.. function:: void fmpq_height(fmpz_t height, const fmpq_t x)

    Sets ``height`` to the height of `x`, defined as the larger of
    the absolute values of the numerator and denominator of `x`.

.. function:: mp_bitcnt_t fmpq_height_bits(const fmpq_t x)

    Returns the number of bits in the height of `x`.


Conversion
--------------------------------------------------------------------------------


.. function:: void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)

    Sets ``res`` to the canonical form of the fraction ``p / q``.
    This is equivalent to assigning the numerator and denominator
    separately and calling ``fmpq_canonicalise``.

.. function:: void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)

    Sets ``a``, ``b`` to the numerator and denominator of ``c``
    respectively.

.. function:: void fmpq_set_si(fmpq_t res, slong p, ulong q)

    Sets ``res`` to the canonical form of the fraction ``p / q``.

.. function:: void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q)

    Sets ``(rnum, rden)`` to the canonical form of the fraction
    ``p / q``. ``rnum`` and ``rden`` may not be aliased.

.. function:: void fmpq_set_mpq(fmpq_t dest, const mpq_t src)

    Sets the value of ``dest`` to that of the ``mpq_t`` variable
    ``src``.

.. function:: void fmpq_set_str(fmpq_t dest, const char * s, int base)

    Sets the value of ``dest`` to the value represented in the string
    ``s`` in base ``base``.

    Returns 0 if no error occurrs. Otherwise returns -1 and ``dest`` is
    set to zero.

.. function:: void fmpq_init_set_mpz_frac_readonly(fmpq_t z, const mpz_t p, const mpz_t q)

    Assuming ``z`` is an ``fmpz_t`` which will not be cleaned up,
    this temporarily copies ``p`` and ``q`` into the numerator and
    denominator of ``z`` for read only operations only. The user must not
    run ``fmpq_clear`` on ``z``.

.. function:: void fmpq_get_mpq(mpq_t dest, const fmpq_t src)

    Sets the value of ``dest``

.. function:: int fmpq_get_mpfr(mpfr_t dest, const fmpq_t src, mpfr_rnd_t rnd)

    Sets the MPFR variable ``dest`` to the value of ``src``,
    rounded to the nearest representable binary floating-point value
    in direction ``rnd``. Returns the sign of the rounding,
    according to MPFR conventions.

.. function:: char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)

.. function:: char * fmpq_get_str(char * str, int b, const fmpq_t x)

    Prints the string representation of `x` in base `b \in [2, 36]` 
    to a suitable buffer.

    If ``str`` is not ``NULL``, this is used as the buffer and 
    also the return value.  If ``str`` is ``NULL``, allocates 
    sufficient space and returns a pointer to the string.

.. function:: void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)

    Sets the uninitialised ``mpq_t`` `z` to the value of the 
    readonly ``fmpq_t`` `f`.

    Note that it is assumed that `f` does not change during 
    the lifetime of `z`.

    The rational `z` has to be cleared by a call to 
    ``flint_mpq_clear_readonly()``.

    The suggested use of the two functions is as follows::

        fmpq_t f;
        ...
        {
            mpq_t z;

            flint_mpq_init_set_readonly(z, f);
            foo(..., z);
            flint_mpq_clear_readonly(z);
        }

    This provides a convenient function for user code, only 
    requiring to work with the types ``fmpq_t`` and ``mpq_t``.

.. function:: void flint_mpq_clear_readonly(mpq_t z)

    Clears the readonly ``mpq_t`` `z`.

.. function:: void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)

    Sets the uninitialised ``fmpq_t`` `f` to a readonly 
    version of the rational `z`.

    Note that the value of `z` is assumed to remain constant 
    throughout the lifetime of `f`.

    The ``fmpq_t`` `f` has to be cleared by calling the 
    function ``fmpq_clear_readonly()``.

    The suggested use of the two functions is as follows::

        mpq_t z;
        ...
        {
            fmpq_t f;

            fmpq_init_set_readonly(f, z);
            foo(..., f);
            fmpq_clear_readonly(f);
        }

.. function:: void fmpq_clear_readonly(fmpq_t f)

    Clears the readonly ``fmpq_t`` `f`.


Input and output
--------------------------------------------------------------------------------


.. function:: int fmpq_fprint(FILE * file, const fmpq_t x)

    Prints ``x`` as a fraction to the stream ``file``.   
    The numerator and denominator are printed verbatim as integers, 
    with a forward slash (/) printed in between.

    In case of success, returns a positive number. In case of failure,
    returns a non-positive number.

.. function:: int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)

    Does the same thing as ``fmpq_fprint``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables. 

    In case of success, returns a positive number. In case of failure,
    returns a non-positive number.

.. function:: int fmpq_print(const fmpq_t x)

    Prints ``x`` as a fraction. The numerator and denominator are
    printed verbatim as integers, with a forward slash (/) printed in
    between.

    In case of success, returns a positive number. In case of failure,
    returns a non-positive number.

.. function:: int _fmpq_print(const fmpz_t num, const fmpz_t den)

    Does the same thing as ``fmpq_print``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables. 

    In case of success, returns a positive number. In case of failure,
    returns a non-positive number.


Random number generation
--------------------------------------------------------------------------------


.. function:: void fmpq_randtest(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)

    Sets ``res`` to a random value, with numerator and denominator
    having up to ``bits`` bits. The fraction will be in canonical
    form. This function has an increased probability of generating
    special values which are likely to trigger corner cases.

.. function:: void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, mp_bitcnt_t bits)

    Does the same thing as ``fmpq_randtest``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    of ``num`` and ``den`` is not allowed.

.. function:: void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)

    As per ``fmpq_randtest``, but the result will not be `0`. 
    If ``bits`` is set to `0`, an exception will result.

.. function:: void fmpq_randbits(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)

    Sets ``res`` to a random value, with numerator and denominator
    both having exactly ``bits`` bits before canonicalisation,
    and then puts ``res`` in canonical form. Note that as a result
    of the canonicalisation, the resulting numerator and denominator can
    be slightly smaller than ``bits`` bits.

.. function:: void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, mp_bitcnt_t bits)

    Does the same thing as ``fmpq_randbits``, but for numerator
    and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    of ``num`` and ``den`` is not allowed.



Arithmetic
--------------------------------------------------------------------------------



.. function:: void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

.. function:: void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

.. function:: void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

.. function:: void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

    Sets ``res`` respectively to ``op1 + op2``, ``op1 - op2``,
    ``op1 * op2``, or ``op1 / op2``. Assumes that the inputs
    are in canonical form, and produces output in canonical form.
    Division by zero results in an error.
    Aliasing between any combination of the variables is allowed.

.. function:: void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

.. function:: void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

.. function:: void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

.. function:: void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

    Sets ``(rnum, rden)`` to the canonical form of the sum,
    difference, product or quotient respectively of the fractions
    represented by ``(op1num, op1den)`` and ``(op2num, op2den)``.
    Aliasing between any combination of the variables is allowed,
    whilst no numerator is aliased with a denominator.

.. function:: void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)

.. function:: void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    
.. function:: void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)

.. function:: void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)

    Sets ``(rnum, rden)`` to the canonical form of the sum or difference
    respectively of the fractions represented by ``(p, q)`` and
    ``(r, 1)``. Numerators may not be aliased with denominators.

.. function:: void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)

.. function:: void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c)

.. function:: void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);

.. function:: void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);

   Sets ``res`` to the sum or difference respectively, of the fraction 
   ``op1`` and the integer `c`.

.. function:: void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

.. function:: void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

    Sets ``res`` to ``res + op1 * op2`` or ``res - op1 * op2``
    respectively, placing the result in canonical form. Aliasing
    between any combination of the variables is allowed.

.. function:: void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

.. function:: void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)

    Sets ``(rnum, rden)`` to the canonical form of the fraction
    ``(rnum, rden)`` + ``(op1num, op1den)`` * ``(op2num, op2den)`` or
    ``(rnum, rden)`` - ``(op1num, op1den)`` * ``(op2num, op2den)``
    respectively. Aliasing between any combination of the variables is allowed,
    whilst no numerator is aliased with a denominator.

.. function:: void fmpq_inv(fmpq_t dest, const fmpq_t src)

    Sets ``dest`` to ``1 / src``. The result is placed in canonical
    form, assuming that ``src`` is already in canonical form.

.. function:: void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, const fmpz_t opnum, const fmpz_t opden, slong e);

.. function:: void fmpq_pow_si(fmpq_t res, const fmpq_t op, slong e);

    Sets ``res`` to ``op`` raised to the power~`e`, where~`e` 
    is a ``slong``.  If `e` is `0` and ``op`` is `0`, then 
    ``res`` will be set to `1`.

.. function:: void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)

    Sets ``res`` to the product of the rational number ``op`` 
    and the integer ``x``.

.. function:: void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)

    Sets ``res`` to the quotient of the rational number ``op`` 
    and the integer ``x``.

.. function:: void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, mp_bitcnt_t exp)

    Sets ``res`` to ``x`` multiplied by ``2^exp``.

.. function:: void fmpq_div_2exp(fmpq_t res, const fmpq_t x, mp_bitcnt_t exp)

    Sets ``res`` to ``x`` divided by ``2^exp``.

.. function:: _fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)

    Set ``(rnum, rden)`` to the gcd of ``(p, q)`` and ``(r, s)``
    which we define to be the canonicalisation of gcd`(ps, qr)/(qs)`. (This is
    apparently Euclid's original definition and is stable under scaling of
    numerator and denominator. It also agrees with the gcd on the integers.
    Note that it does not agree with gcd as defined in ``fmpq_poly``.)
    This definition agrees with the result as output by Sage and Pari/GP.

.. function:: fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)

    Set ``res`` to the gcd of ``op1`` and ``op2``. See the low
    level function ``_fmpq_gcd`` for our definition of gcd.


Modular reduction and rational reconstruction
--------------------------------------------------------------------------------


.. function:: int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)

.. function:: int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)

    Sets the integer ``res`` to the residue `a` of
    `x = n/d` = ``(num, den)`` modulo the positive integer `m` = ``mod``,
    defined as the `0 \le a < m` satisfying `n \equiv a d \pmod m`.
    If such an `a` exists, 1 will be returned, otherwise 0 will
    be returned.

.. function:: int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)

.. function:: int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)

    Reconstructs a rational number from its residue `a` modulo `m`.

    Given a modulus `m > 1`, a residue `0 \le a < m`, and positive `N, D`
    satisfying `2ND < m`, this function attempts to find a fraction `n/d` with
    `0 \le |n| \le N` and `0 < d \le D` such that `\gcd(n,d) = 1` and
    `n \equiv ad \pmod m`. If a solution exists, then it is also unique.
    The function returns 1 if successful, and 0 to indicate that no solution
    exists.

.. function:: int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)

.. function:: int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)

    Reconstructs a rational number from its residue `a` modulo `m`,
    returning 1 if successful and 0 if no solution exists.
    Uses the balanced bounds `N = D = \lfloor\sqrt{m/2}\rfloor`.



Rational enumeration
--------------------------------------------------------------------------------


.. function:: void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

.. function:: void fmpq_next_minimal(fmpq_t res, const fmpq_t x)

    Given `x` which is assumed to be nonnegative and in canonical form, sets
    ``res`` to the next rational number in the sequence obtained by
    enumerating all positive denominators `q`, for each `q` enumerating
    the numerators `1 \le p < q` in order and generating both `p/q` and `q/p`,
    but skipping all `\gcd(p,q) \ne 1`. Starting with zero, this generates
    every nonnegative rational number once and only once, with the first
    few entries being:

    ``0, 1, 1/2, 2, 1/3, 3, 2/3, 3/2, 1/4, 4, 3/4, 4/3, 1/5, 5, 2/5, \ldots.``

    This enumeration produces the rational numbers in order of
    minimal height. It has the disadvantage of being somewhat slower to
    compute than the Calkin-Wilf enumeration.

.. function:: void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

.. function:: void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)

    Given a signed rational number `x` assumed to be in canonical form, sets
    ``res`` to the next element in the minimal-height sequence
    generated by ``fmpq_next_minimal`` but with negative numbers
    interleaved:

    ``0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots.``

    Starting with zero, this generates every rational number once
    and only once, in order of minimal height.

.. function:: void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

.. function:: void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)

    Given `x` which is assumed to be nonnegative and in canonical form, sets
    ``res`` to the next number in the breadth-first traversal of the
    Calkin-Wilf tree. Starting with zero, this generates every nonnegative
    rational number once and only once, with the first few entries being:

    ``0, 1, 1/2, 2, 1/3, 3/2, 2/3, 3, 1/4, 4/3, 3/5, 5/2, 2/5, \ldots.``

    Despite the appearance of the initial entries, the Calkin-Wilf
    enumeration does not produce the rational numbers in order of height:
    some small fractions will appear late in the sequence. This order
    has the advantage of being faster to produce than the minimal-height
    order.

.. function:: void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)

.. function:: void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)

    Given a signed rational number `x` assumed to be in canonical form, sets
    ``res`` to the next element in the Calkin-Wilf sequence with
    negative numbers interleaved:

    ``0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots.``

    Starting with zero, this generates every rational number once
    and only once, but not in order of minimal height.


Continued fractions
--------------------------------------------------------------------------------


.. function:: slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)

    Generates up to `n` terms of the (simple) continued fraction expansion
    of `x`, writing the coefficients to the vector `c` and the remainder `r`
    to the ``rem`` variable. The return value is the number `k` of
    generated terms. The output satisfies

    .. math ::

        x = c_0 + \cfrac{1}{c_1 + \cfrac{1}{c_2 +
            \cfrac{1}{ \ddots + \cfrac{1}{c_{k-1} + r }}}}

    If `r` is zero, the continued fraction expansion is complete.
    If `r` is nonzero, `1/r` can be passed back as input to generate
    `c_k, c_{k+1}, \ldots`. Calls to ``fmpq_get_cfrac`` can therefore
    be chained to generate the continued fraction incrementally,
    extracting any desired number of coefficients at a time.

    In general, a rational number has exactly two continued fraction
    expansions. By convention, we generate the shorter one. The longer
    expansion can be obtained by replacing the last coefficient
    `a_{k-1}` by the pair of coefficients `a_{k-1} - 1, 1`.

    As a special case, the continued fraction expansion of zero consists
    of a single zero (and not the empty sequence).

    This function implements a simple algorithm, performing repeated
    divisions. The running time is quadratic.

.. function:: void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n)

    Sets `x` to the value of the continued fraction

    .. math ::

        x = c_0 + \cfrac{1}{c_1 + \cfrac{1}{c_2 +
            \cfrac{1}{ \ddots + \cfrac{1}{c_{n-1}}}}}

    where all `c_i` except `c_0` should be nonnegative.
    It is assumed that `n > 0`.

    For large `n`, this function implements a subquadratic algorithm.
    The convergents are given by a chain product of 2 by 2 matrices.
    This product is split in half recursively to balance the size
    of the coefficients.

.. function:: slong fmpq_cfrac_bound(const fmpq_t x)

    Returns an upper bound for the number of terms in the continued
    fraction expansion of `x`. The computed bound is not necessarily sharp.

    We use the fact that the smallest denominator
    that can give a continued fraction of length `n` is the Fibonacci
    number `F_{n+1}`.


Special functions
--------------------------------------------------------------------------------


.. function:: void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n)

.. function:: void fmpq_harmonic_ui(fmpq_t x, ulong n)

    Computes the harmonic number `H_n = 1 + 1/2 + 1/3 + \dotsb + 1/n`.
    Table lookup is used for `H_n` whose numerator and denominator 
    fit in single limb. For larger `n`, a divide and conquer strategy is used.


Dedekind sums
--------------------------------------------------------------------------------

Most of the definitions and relations used in the following section
are given by Apostol \cite{Apostol1997}. The Dedekind sum `s(h,k)` is
defined for all integers `h` and `k` as

.. math ::

    s(h,k) = \sum_{i=1}^{k-1} \left(\left(\frac{i}{k}\right)\right)
    \left(\left(\frac{hi}{k}\right)\right)

where 

.. math ::

    ((x))=\begin{cases}
    x-\lfloor x\rfloor-1/2 &\mbox{if }
    x\in\mathbf{Q}\setminus\mathbf{Z}\\
    0 &\mbox{if }x\in\mathbf{Z}.
    \end{cases}

If `0 < h < k` and `(h,k) = 1`, this reduces to

.. math ::

    s(h,k) = \sum_{i=1}^{k-1} \frac{i}{k}
        \left(\frac{hi}{k}-\left\lfloor\frac{hi}{k}\right\rfloor
        -\frac{1}{2}\right).

The main formula for evaluating the series above is the following.
Letting `r_0 = k`, `r_1 = h`, `r_2, r_3, \ldots, r_n, r_{n+1} = 1`
be the remainder sequence in the Euclidean algorithm for
computing GCD of `h` and `k`, 
`s(h,k) = \frac{1-(-1)^n}{8} - \frac{1}{12} \sum_{i=1}^{n+1}
(-1)^i \left(\frac{1+r_i^2+r_{i-1}^2}{r_i r_{i-1}}\right).`
Writing `s(h,k) = p/q`, some useful properties employed are
`|s| < k / 12`, `q | 6k` and `2|p| < k^2`.

.. function:: void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)

    Computes `s(h,k)` for arbitrary `h` and `k` using a straightforward
    implementation of the defining sum using ``fmpz`` arithmetic.
    This function is slow except for very small `k` and is mainly
    intended to be used for testing purposes.

.. function:: double fmpq_dedekind_sum_coprime_d(double h, double k)

    Returns an approximation of `s(h,k)` computed by evaluating the
    remainder sequence sum using double-precision arithmetic.
    Assumes that `0 < h < k` and `(h,k) = 1`, and that `h`, `k` and
    their remainders can be represented exactly as doubles, e.g.
    `k < 2^{53}`.

    We give a rough error analysis with IEEE double precision arithmetic,
    assuming `2 k^2 < 2^{53}`. By assumption, the terms in the sum evaluate
    exactly apart from the division. Thus each term is bounded in magnitude
    by `2k` and its absolute error is bounded by `k 2^{-52}`.
    By worst-case analysis of the Euclidean algorithm, we also know that
    no more than 40 terms will be added.

    It follows that the absolute error is at most `C k 2^{-53}` for
    some constant `C`. If we multiply the output by `6 k` in order
    to obtain an integer numerator, the order of magnitude of the error
    is around `6 C k^2 2^{-53}`, so rounding to the nearest integer gives
    a correct numerator whenever `k < 2^{26-d}` for some small number of
    guard bits `d`. A computation has shown that `d = 5` is sufficient,
    i.e. this function can be used for exact computation when
    `k < 2^{21} \approx 2 \times 10^6`. This bound can likely be improved.

.. function:: void fmpq_dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k)

    Computes `s(h,k)` for `h` and `k` satisfying `0 \le h \le k` and
    `(h,k) = 1`. This function effectively evaluates the remainder
    sequence sum using ``fmpz`` arithmetic, without optimising for
    any special cases. To avoid rational arithmetic, we use
    the integer algorithm of Knuth \cite{Knuth1977}.

.. function:: void fmpq_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k)

    Computes `s(h,k)` for `h` and `k` satisfying `0 \le h \le k`
    and `(h,k) = 1`.

    This function calls ``fmpq_dedekind_sum_coprime_d`` if `k` is small
    enough for a double-precision estimate of the sum to yield a correct
    numerator upon multiplication by `6k` and rounding to the nearest integer.
    Otherwise, it calls ``fmpq_dedekind_sum_coprime_large``.

.. function:: void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)

    Computes `s(h,k)` for arbitrary `h` and `k`. If the caller
    can guarantee `0 < h < k` and `(h,k) = 1` ahead of time, it is always
    cheaper to call ``fmpq_dedekind_sum_coprime``.

    This function uses the following identities to reduce the general
    case to the situation where `0 < h < k` and `(h,k) = 1`:
    If `k \le 2` or `h = 0`, `s(h,k) = 0`.
    If `h < 0`, `s(h,k) = -s(-h,k)`.
    For any `q > 0`, `s(qh,qk) = s(h,k)`.
    If `0 < k < h` and `(h,k) = 1`,
    `s(h,k) = (1+h(h-3k)+k^2) / (12hk) - t(k,h).`

