.. _gr-mpoly:

**gr_mpoly.h** -- sparse multivariate polynomials over generic rings
===============================================================================

This module implements multivariate polynomials
with ``gr`` coefficients.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_mpoly_struct
          gr_mpoly_t

    Represents a multivariate polynomial
    `f \in R[X_1,\ldots,X_n]` as an array of coefficients
    in a generic ring *R* together with an array of packed exponents.
    The two arrays always have the same length.

    A :type:`gr_mpoly_t` is always normalised by removing zero
    coefficients.
    For rings without decidable equality (e.g. rings with inexact
    representation), only coefficients that are provably zero will be
    removed, and there can thus be spurious zeros in the
    internal representation. For example, with ball coefficients
    one can have the polynomial `3 x y^2 + [\pm 0.01] x y`.
    Over rings with this issue, the represented lengths or degrees are
    thus upper bounds, and methods that depend on knowing the exact
    term structure of a polynomial will return ``GR_UNABLE``
    when encountering such input.

    A ``gr_mpoly_t`` is defined as an array of length one of type
    ``gr_mpoly_struct``, permitting a ``gr_mpoly_t`` to
    be passed by reference.

.. type:: gr_mpoly_ctx_struct
          gr_mpoly_ctx_t

    Context object representing a multivariate polynomial ring
    `R[X_1,\ldots,X_n]`.
    This subtypes :type:`gr_ctx_t`, allowing
    generic ``gr`` and ``gr_ctx`` methods to be used interchangeably
    with ``gr_mpoly`` and ``gr_mpoly_ctx`` methods. For example,
    :func:`gr_add` with :type:`gr_mpoly_t` and :type:`gr_mpoly_ctx_t`
    arguments is equivalent to :func:`gr_mpoly_add`.
    A context object contains the following data:

    * A pointer *mctx* to a :type:`gr_ctx_t` representing the coefficient type *R*.
      The coefficient context object is not considered owned by
      the ``gr_mpoly_ctx`` and the user must ensure
      that it stays alive as long as the ``gr_mpoly_ctx`` is alive.

    * A pointer *cctx* to a :type:`mpoly_ctx_t` defining the number of
      variables and term ordering. This object is considered
      owned by and will automatically be initialized and cleared
      along with the ``gr_mpoly_ctx``.

    * An optional pointer *vars* to an array of strings
      specifying names of the generators `X_1, \ldots, X_n`.
      This can be set with
      :func:`gr_mpoly_ctx_set_gen_names`. By default, *vars* will be
      initialized to ``NULL`` in which case some default names
      are used.
      Names are used for printing and parsing from strings
      with :func:`gr_set_str` and in some cases for coercions
      between different rings.

.. macro:: GR_MPOLY_MCTX(ctx)

    Access the mpoly context object *mctx*.

.. macro:: GR_MPOLY_CCTX(ctx)

    Access the coefficient context object *cctx*.

.. macro:: GR_MPOLY_VARS(ctx)

    Access the array of variable names *vars*.

.. macro:: GR_MPOLY_NVARS(ctx)

    Access the number of variables of this context object.

Context object methods
-------------------------------------------------------------------------------

.. function:: void gr_mpoly_ctx_init(gr_mpoly_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)

    Initializes ``ctx`` to represent a polynomial ring with
    coefficients in ``base_ring``, with ``nvars`` variables
    and term ordering ``ord``.

.. function:: void gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx)

    Clears the context object ``ctx``.

.. function:: void gr_mpoly_ctx_init_rand(gr_mpoly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring, slong max_nvars)

    Initializes ``ctx`` with a random number of variables
    up to ``max_nvars`` inclusive and with a random term ordering.

The following methods implement parts of the standard interface
for ``gr`` context objects.

.. function:: int gr_mpoly_ctx_set_gen_names(gr_mpoly_ctx_t ctx, const char ** s)

    Sets the names of the generators to the strings in ``s``.

.. function:: int gr_mpoly_ctx_write(gr_stream_t out, gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_ring(gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_zero_ring(gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_commutative_ring(gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_integral_domain(gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_field(gr_mpoly_ctx_t ctx)
              truth_t gr_mpoly_ctx_is_threadsafe(gr_mpoly_ctx_t ctx)

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_mpoly_init(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Initializes and sets *A* to the zero polynomial.

.. function:: void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)
              void gr_mpoly_init2(gr_mpoly_t A, slong alloc, gr_mpoly_ctx_t ctx)

    Initializes *A* with space allocated for the given number
    of coefficients and exponents with the given number of bits.

.. function:: void gr_mpoly_clear(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Clears *A*, freeing all allocated data.

Vectors of polynomials
-------------------------------------------------------------------------------

.. type:: gr_mpoly_vec_struct
          gr_mpoly_vec_t

    A resizable vector of :type:`gr_mpoly_t` elements. It is safe to
    cast a :type:`gr_mpoly_vec_t` to a :type:`gr_vec_t` (with :type:`gr_mpoly_t`
    elements of the correct type) and vice versa.

.. function:: void gr_mpoly_vec_init(gr_mpoly_vec_t vec, slong len, gr_ctx_t ctx)

    Initializes *vec* to a vector of length *len* with all entries
    set to the zero polynomial of type *ctx*. The length must be nonnegative.

.. function:: void gr_mpoly_vec_clear(gr_mpoly_vec_t vec, gr_ctx_t ctx)

    Clears the vector *vec*, freeing all allocated memory.

.. function:: gr_mpoly_struct * gr_mpoly_vec_entry_ptr(gr_mpoly_vec_t vec, slong i, gr_ctx_t ctx)
              const gr_mpoly_struct * gr_mpoly_vec_entry_srcptr(const gr_mpoly_vec_t vec, slong i, gr_ctx_t ctx)

    Returns a pointer to the *i*-th polynomial in the vector, indexed from zero.
    The index must be in bounds.

.. function:: slong gr_mpoly_vec_length(const gr_mpoly_vec_t vec, gr_ctx_t ctx)

    Returns the length of the vector.

.. function:: void gr_mpoly_vec_fit_length(gr_mpoly_vec_t vec, slong len, gr_ctx_t ctx)

    Allocates space for at least *len* elements. This does not change
    the length of the vector.

.. function:: void gr_mpoly_vec_set_length(gr_mpoly_vec_t vec, slong len, gr_ctx_t ctx)

    Resizes the vector to length *len*, which must be nonnegative.
    The vector will be extended with zero polynomials if necessary.

.. function:: int gr_mpoly_vec_set(gr_mpoly_vec_t res, const gr_mpoly_vec_t src, gr_ctx_t ctx)

    Sets *res* to a copy of *src*.

.. function:: int gr_mpoly_vec_append(gr_mpoly_vec_t vec, const gr_mpoly_t f, gr_ctx_t ctx)

    Appends the polynomial *f* to the end of the vector.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void _gr_mpoly_normalise(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Removes provably zero coefficients from ``A`` and updates the length.
    If all coefficients are zero, the length is set to zero.  This function
    is mainly used internally, as all functions guarantee normalisation.

.. function:: void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Swaps *A* and *B* efficiently.

.. function:: void gr_mpoly_set_shallow(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Sets *A* to a shallow copy of *B* (unsafe).

.. function:: int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Sets *A* to *B*.

.. function:: int gr_mpoly_zero(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Sets *A* to the zero polynomial.

.. function:: truth_t gr_mpoly_is_zero(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Returns whether *A* is the zero polynomial.

.. function:: slong gr_mpoly_length(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Returns the number of terms in *A*.

Generators
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_gen(gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx)

    Sets *A* to the generator with index *var* (indexed from zero).

.. function:: truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx)

    Returns whether *A* is the generator with index *var* (indexed from zero).

.. function:: int gr_mpoly_gens(gr_vec_t res, gr_mpoly_ctx_t ctx)

    Sets the vector *res* to a list of the generators `X_1, \ldots, X_n`.

.. function:: int gr_mpoly_gens_recursive(gr_vec_t vec, gr_mpoly_ctx_t ctx)

    Sets the vector *res* to a list of the recursive generators of `R`
    (as constant elements of `R[X_1, \ldots, X_n]`)
    followed by the generators `X_1, \ldots, X_n`.

Conversions
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_set_scalar(gr_mpoly_t A, gr_srcptr c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_ui(gr_mpoly_t A, ulong c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_si(gr_mpoly_t A, slong c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_fmpz(gr_mpoly_t A, const fmpz_t c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_fmpq(gr_mpoly_t A, const fmpq_t c, gr_mpoly_ctx_t ctx)

    Sets *A* to the given scalar *c*.

.. function:: int gr_mpoly_set_other(gr_mpoly_t res, gr_srcptr A, gr_ctx_t A_ctx, gr_mpoly_ctx_t ctx)

    Sets *res* to *A* (an element of *A_ctx*) converted to the
    multivariate polynomial ring *ctx*.

    If *A_ctx* is a multivariate polynomial ring, this attempts to
    coerce the coefficients and translate the generators.
    If both rings have named generators, we find all used
    generators in *A* and match them to generators with the same names
    in *ctx*. If both rings have the same number of unnamed generators
    and the same term ordering, we perform a direct conversion.
    Other cases are not currently supported.
 
    Otherwise, we attempt to interpret *A* as a scalar.

    Currently, absorbing generators from nested rings is not supported,
    e.g. converting between `R[x,y][s,t]` and `R[x,y,s,t]` is likely to fail.


Comparisons
-------------------------------------------------------------------------------

.. function:: truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Returns whether *A* and *B* are equal.

Random generation
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, gr_mpoly_ctx_t ctx)

    Sets *A* to a random polynomial with up to *length* terms
    and up to *exp_bits* bits in the exponents.

Input and output
-------------------------------------------------------------------------------

Note: :func:`gr_set_str` can be used for parsing.

.. function:: int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
              int gr_mpoly_print_pretty(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

    Prints *A*, using the variable names stored in the context.

Coefficient and exponent access
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx)

    Sets *c* to the coefficient in *A* with exponents *exp*.

.. function:: int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, gr_mpoly_ctx_t ctx)

.. function:: int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, gr_mpoly_ctx_t ctx)
              int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, gr_mpoly_ctx_t ctx)

    Sets the coefficient with exponents *exp* in *A* to the scalar *c*
    which must be an element of or coercible to the coefficient ring.

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Sets *A* to the negation of *B*.

.. function:: int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)

    Sets *A* to the difference of *B* and *C*.

.. function:: int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)

    Sets *A* to the difference of *B* and *C*.

.. function:: int gr_mpoly_mul(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_heap(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_heap_threaded(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)

    Sets *A* to the product of *B* and *C*.
    The *monomial* version assumes that *C* is a monomial.
    The *heap* and *heap_threaded* versions implement Johnson's heap-based
    algorithms.

.. function:: int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, gr_mpoly_ctx_t ctx)
              int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, gr_mpoly_ctx_t ctx)

    Sets *A* to *B* multiplied by the scalar *c* which must be
    an element of or coercible to the coefficient ring.

.. function:: int gr_mpoly_sqr_monomial(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
              int gr_mpoly_sqr_commutative_heap(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
              int gr_mpoly_sqr_commutative_heap_threaded(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
              int gr_mpoly_sqr(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)

    Sets *poly1* to the square of *poly2*. The *monomial* version assumes
    that *poly2* is a monomial. The *commutative_heap* and
    *commutative_heap_threaded* versions assume a commutative coefficient ring.

Division
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_inv(gr_mpoly_t res, const gr_mpoly_t src, gr_mpoly_ctx_t ctx)

    If *src* has a multiplicative inverse, sets *res* to the inverse
    and returns ``GR_SUCCESS``. If *src* is not invertible,
    returns ``GR_DOMAIN``.
    If invertibility cannot be decided, returns ``GR_UNABLE``.

.. function:: int gr_mpoly_divides_heap(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_divides(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Standard checked, exact division.
    If *A* is divisible by *B*, sets *Q* to the exact quotient `Q = A / B` and returns ``GR_SUCCESS``.
    If *A* is not divisible by *B*, sets *Q* to zero and returns ``GR_DOMAIN``.
    If divisibility cannot be decided, returns ``GR_UNABLE``.

    The *heap* version implements the Monagan-Pearce heap algorithm.
    The default version currently only delegates to the heap algorithm
    except for trivial special cases.

.. function:: int gr_mpoly_divrem_heap(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_divrem(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Standard multivariate division with remainder, typically used
    when ``ctx`` is a field.
    Sets *Q* and *R* to the quotient and remainder of *A* divided by *B*,
    satisfying `A = Q B + R` and where no monomial in *R* is divisible
    by the leading monomial of *B*.

    If ``ctx`` is not a field, this works if the leading coefficient in *B*
    is a unit and more generally if any division by the leading coefficient in *B*
    succeeds in the execution of the division algorithm. If any non-exact
    division is encountered, returns ``GR_DOMAIN``.

.. function:: int gr_mpoly_div(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

    Like :func:`gr_mpoly_divrem`, but discarding the remainder.

.. function:: int gr_mpoly_divrem_ideal(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)

    Like :func:`gr_mpoly_divrem`, but takes a vector *B* of divisor polynomials
    and outputs a vector *Q* of quotient polynomials such that
    `A = \sum_i Q_i B_i + R`, where no monomial in *R* is divisible by
    any leading monomial in *B*.

.. function:: int gr_mpoly_div_weak(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_divrem_weak(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_divrem_ideal_weak(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)

    Like the ``divrem`` functions, producing quotient and remainder
    satisfying `A = Q B + R` (or `A = \sum_i Q_i B_i + R` in the ideal case),
    except that if a nonexact coefficient division is encountered,
    instead of halting and returning ``GR_DOMAIN``,
    we divide the coefficients with remainder and accumulate the
    remainder. This produces a polynomial remainder which is partially
    reduced, but where monomials in *R* may still be divisible (ignoring
    the coefficient) by the leading monomial of *B*.

.. function:: int gr_mpoly_quasidivrem_heap(gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_quasidiv(gr_ptr scale, gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_quasidivrem(gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
              int gr_mpoly_quasidivrem_ideal(gr_ptr scale, gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)

    Like the ``divrem`` functions, but compute a scalar ``scale`` (as an
    element of the scalar ring of ``ctx``) such that
    `scale \cdot A = Q B + R` (or `scale \cdot A = \sum_i Q_i B_i + R` in the ideal case).
    This is typically useful to perform division over the fraction field when ``ctx`` is an
    integral domain.

Derivative and integral
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_derivative(gr_mpoly_t A, const gr_mpoly_t B, slong var, gr_mpoly_ctx_t ctx)

    Set *A* to the derivative of *B* with respect to the variable of index *var*.

.. function:: int gr_mpoly_integral(gr_mpoly_t A, const gr_mpoly_t B, slong var, gr_mpoly_ctx_t ctx)

    Set *A* to the integral of *B* with respect to the variable of index *var*.

Other operations
-------------------------------------------------------------------------------

.. function:: int gr_mpoly_canonical_associate(gr_mpoly_t res, gr_mpoly_t u, const gr_mpoly_t src, gr_mpoly_ctx_t ctx)


Container operations
-------------------------------------------------------------------------------

Mostly intended for internal use.

.. function:: void _gr_mpoly_fit_length(gr_ptr * coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length, gr_mpoly_ctx_t ctx)

.. function:: void gr_mpoly_fit_length(gr_mpoly_t A, slong len, gr_mpoly_ctx_t ctx)

    Ensures that *A* has space for *len* coefficients and exponents.

.. function:: void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)

.. function:: void gr_mpoly_fit_length_fit_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)

.. function:: void gr_mpoly_fit_length_reset_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)

.. function:: void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen, gr_mpoly_ctx_t ctx)

.. function:: void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx)

.. function:: int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx)

.. function:: void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx)

.. function:: int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx)

.. function:: void gr_mpoly_sort_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

.. function:: int gr_mpoly_combine_like_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

.. function:: truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

.. function:: void gr_mpoly_assert_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)



.. raw:: latex

    \newpage
