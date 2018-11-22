.. _fmpz-mpoly:

**fmpz_mpoly.h** -- multivariate polynomials over the integers
===============================================================================

Description.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mpoly_struct

.. type:: fmpz_mpoly_t

    Description.

.. type:: fmpz_mpoly_ctx_struct

.. type:: fmpz_mpoly_ctx_t

    Description.

.. type:: fmpz_mpoly_univar_struct

.. type:: fmpz_mpoly_univar_t

    Description.


Context object
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

    Initialise a context object for a polynomial ring with the given number of
    variables and the given ordering. The possibilities for the ordering are
    ``ORD_LEX``, ``ORD_REVLEX``, ``ORD_DEGLEX`` and
    ``ORD_DEGREVLEX``.

.. function:: slong fmpz_mpoly_ctx_nvars(fmpz_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)

    Release up any space allocated by an ``fmpz_mpoly_ctx_t``.


Memory management
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_init(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Initialise an ``fmpz_mpoly_t`` for use, given an initialised context
    object. Its value is set to 0.
    By default 8 bits are allocated for the exponent widths.

.. function:: void fmpz_mpoly_init2(fmpz_mpoly_t poly, slong alloc, const fmpz_mpoly_ctx_t ctx)

    Initialise an ``fmpz_mpoly_t`` for use, with space for at least
    ``alloc`` terms, given an initialised context. Its value is set to 0.
    By default 8 bits are allocated for the exponent widths.

.. function:: void fmpz_mpoly_realloc(fmpz_mpoly_t poly, slong len, const fmpz_mpoly_ctx_t ctx)

    Reallocate an ``fmpz_mpoly_t`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than
    ``len``.

.. function:: void _fmpz_mpoly_fit_length(fmpz ** poly, ulong ** exps, slong * alloc, slong len, slong N)

    Reallocate a low level ``fmpz_mpoly`` to have space for at least
    ``len`` terms. No truncation is performed if ``len`` is less than
    the currently allocated number of terms; the allocated space can only grow.
    Assumes exponent vectors each consist of `N` words.

.. function:: void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, const fmpz_mpoly_ctx_t ctx)

    Reallocate a low level ``fmpz_mpoly`` to have space for at least
    ``len`` terms. No truncation is performed if ``len`` is less than
    the currently allocated number of terms; the allocated space can only grow.

.. function:: void _fmpz_mpoly_set_length(fmpz_mpoly_t poly, slong newlen, const fmpz_mpoly_ctx_t ctx)

    Set the number of terms of the given polynomial to the given length. 
    Assumes the polynomial has at least ``newlen`` allocated and initialised
    terms.

.. function:: void fmpz_mpoly_fit_bits(fmpz_mpoly_t poly, slong bits, const fmpz_mpoly_ctx_t ctx)

    Reallocate the polynomial to have space for exponent fields of the given
    number of bits. This function can increase the number of bits only.

.. function:: void fmpz_mpoly_clear(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Release any space allocated for an ``fmpz_mpoly_t``.


Basic manipulation
----------------------------------------------------------------------


.. function:: int _fmpz_mpoly_fits_small(const fmpz * poly, slong len)

    Return 1 if the array of coefficients of length ``len`` consists
    entirely of values that are small ``fmpz`` values, i.e. of at most
    ``FLINT_BITS - 2`` bits plus a sign bit.

.. function:: slong fmpz_mpoly_max_bits(const fmpz_mpoly_t poly)

    Computes the maximum number of bits `b` required to represent the absolute
    values of the coefficients of ``poly``. If all of the coefficients are
    positive, `b` is returned, otherwise `-b` is returned.



Degrees
----------------------------------------------------------------------


.. function:: int fmpz_mpoly_degrees_fit_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if the degrees of ``A`` with respect to each variable fit into an ``slong``, otherwise return ``0``.

.. function:: void fmpz_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``A`` with respect to each variable.
    If ``A`` is zero, all degrees are set to ``-1``.

.. function:: void fmpz_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

    Either return or set ``deg`` to the degree of ``A`` with respect to the variable of index ``var``.
    If ``A`` is zero, the degree is defined to be ``-1``.

.. function:: int fmpz_mpoly_total_degree_fits_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if the total degree of ``A`` fits into an ``slong``, otherwise return ``0``.

.. function:: void fmpz_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Either return or set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, the total degree is defined to be ``-1``.


Coefficients
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t poly, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``poly2`` is a monomial,
    set `c` to the coefficient of the corresponding monomial in ``poly``.
    This function thows if ``poly2`` is not a monomial.

.. function:: void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t poly, const fmpz_t c, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``poly2`` is a monomial,
    set the coefficient of the corresponding monomial in ``poly`` to `c`.
    This function thows if ``poly2`` is not a monomial.

.. function:: void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t poly, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

    Set `c` to the coefficient of the monomial with exponent vector ``exp``.

.. function:: void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly, const fmpz_t c, fmpz * const * exp, fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.

.. function:: void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t poly, ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.

.. function:: void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t poly, slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.


.. function:: void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t poly, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

    Set `c` to the coefficient of the monomial with exponent vector ``exp``.

.. function:: void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t poly, const fmpz_t c, ulong const * exp, fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.

.. function:: void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t poly, ulong c, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.

.. function:: void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t poly, slong c, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to `c`.


Container operations
----------------------------------------------------------------------

    These functions deal with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.

.. function:: int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is in canonical form. Otherwise, return ``0``.
    To be in canonical form, all of the terms must have nonzero coefficient, and the terms must be sorted from greatest to least.

.. function:: slong fmpz_mpoly_length(const fmpz_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return the number of terms in ``A``.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

.. function:: ulong fmpz_mpoly_get_term_coeff_ui(const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_get_term_coeff_si(const fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)

    Either return or set ``c`` to the coefficient of the term of index ``i``.

.. function:: void fmpz_mpoly_set_term_coeff_fmpz(fmpz_mpoly_t A, slong i, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_term_coeff_ui(fmpz_mpoly_t A, slong i, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_term_coeff_si(fmpz_mpoly_t A, slong i, slong c, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the term of index ``i`` to ``c``.

.. function:: int fmpz_mpoly_term_exp_fits_si(const fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)

.. function:: int fmpz_mpoly_term_exp_fits_ui(const fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if all entries of the exponent vector of the term of index `i`  fit into an ``slong`` (resp. a ``ulong``). Otherwise, return ``0``.

.. function:: void fmpz_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Set ``exp`` to the exponent vector of the term of index ``i``.

.. function:: void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_termexp_ui(fmpz_mpoly_t A, slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    Set the exponent vector of the term of index ``i`` to ``exp``.

.. function:: void fmpz_mpoly_push_term_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_push_term_ui_fmpz(fmpz_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_push_term_si_fmpz(fmpz_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A, ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A, slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    Append a term to ``A`` with coefficient ``c`` and exponent vector ``exp``.
    This function runs in constant average time.

.. function:: void fmpz_mpoly_sort_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Sort the terms of ``A`` into the canonical ordering dictated by the ordering in ``ctx``.
    This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    This function runs in linear time in the bit size of ``A``.

.. function:: void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Combine adjacent like terms in ``A`` and delete terms with coefficient zero.
    If the terms of ``A`` were sorted to begin with, the result will be in canonical form.
    This function runs in linear time in the bit size of ``A``.

.. function:: void fmpz_mpoly_reverse(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the reversal of ``B``.


Set and negate
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
    
    Set ``poly1`` to ``poly2``.

.. function:: void fmpz_mpoly_swap(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Efficiently swap the contents of the two given polynomials. No copying is
    performed; the swap is accomplished by swapping pointers.

.. function:: void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the `i`-th generator (variable),
    where `i = 0` corresponds to the variable with the most significance
    with respect to the ordering. 

.. function:: void fmpz_mpoly_neg(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
    
    Set ``poly1`` to ``-poly2``.


Constants
----------------------------------------------------------------------


.. function:: int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is a constant, else return 0.

.. function:: void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``poly`` is a constant, set `c` to this constant.
    This function throws if ``poly`` is not a constant.

.. function:: void fmpz_mpoly_set_fmpz(fmpz_mpoly_t poly, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the constant `c`.

.. function:: void fmpz_mpoly_set_ui(fmpz_mpoly_t poly, ulong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the constant `c`.

.. function:: void fmpz_mpoly_set_si(fmpz_mpoly_t poly, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the constant `c`.

.. function:: void fmpz_mpoly_zero(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the constant 0.

.. function:: void fmpz_mpoly_one(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Set ``poly`` to the constant 1.


Comparison
----------------------------------------------------------------------


.. function:: int fmpz_mpoly_equal(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly1`` is equal to ``poly2``, else return 0.

.. function:: int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t poly, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpz_mpoly_equal_ui(const fmpz_mpoly_t poly, ulong  c, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpz_mpoly_equal_si(const fmpz_mpoly_t poly, slong  c, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpz_mpoly_is_zero(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant 0, else return 0.

.. function:: int fmpz_mpoly_is_one(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant 1, else return 0.


.. function:: int fmpz_mpoly_is_gen(const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    If `i \ge 0`, return 1 if ``poly`` is equal to the `i`-th generator,
    otherwise return 0. If `i < 0`, return 1 if the polynomial is
    equal to any generator, otherwise return 0.


Basic arithmetic
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_add_fmpz(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpz_mpoly_add_ui(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpz_mpoly_add_si(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpz_mpoly_sub_ui(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpz_mpoly_sub_si(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpz_mpoly_add(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus ``poly3``.

.. function:: void fmpz_mpoly_sub(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus ``poly3``.


Scalar operations
----------------------------------------------------------------------


.. function:: void _fmpz_mpoly_scalar_mul_ui(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, ulong c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` times the
    unsigned integer `c`. The exponents vectors are assumed to each consist of
    `N` words. The ouput polynomial is assumed to have space for ``len``
    terms.

.. function:: void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times the unsigned integer `c`.

.. function:: void _fmpz_mpoly_scalar_mul_si(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, slong c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` times the
    signed integer `c`. The exponents vectors are assumed to each consist of
    `N` words. The ouput polynomial is assumed to have space for ``len``
    terms.

.. function:: void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times the signed integer `c`.

.. function:: void _fmpz_mpoly_scalar_mul_fmpz(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, fmpz_t c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` times the
    multiprecision integer `c`. The exponents vectors are assumed to each
    consist of `N` words. The ouput polynomial is assumed to have space for
    ``len`` terms.

.. function:: void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times the multiprecision integer `c`.

.. function:: void _fmpz_mpoly_scalar_divexact_ui(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, ulong c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` divided by the
    unsigned integer `c`. The exponents vectors are assumed to each consist of
    `N` words. The ouput polynomial is assumed to have space for ``len``
    terms. The division is assumed to be exact.

.. function:: void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by the unsigned integer `c`. The
    division is assumed to be exact.

.. function:: void _fmpz_mpoly_scalar_divexact_si(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, slong c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` divided by the
    signed integer `c`. The exponents vectors are assumed to each consist of
    `N` words. The ouput polynomial is assumed to have space for ``len``
    terms. The division is assumed to be exact.

.. function:: void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by the signed integer `c`. The
    division is assumed to be exact.

.. function:: void _fmpz_mpoly_scalar_divexact_fmpz(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len, slong N, fmpz_t c)

    Set ``(poly1, exps1, len)`` to ``(poly2, exps2, len)`` divided by the
    multiprecision integer `c`. The exponents vectors are assumed to each
    consist of `N` words. The ouput polynomial is assumed to have space for
    ``len`` terms. The division is assumed to be exact.

.. function:: void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by the multiprecision integer `c`.
    The division is assumed to be exact.


Multiplication
----------------------------------------------------------------------


.. function:: slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong N)

    Set ``(poly1, exp1, alloc)`` to ``(poly2, exps2, len2)`` times
    ``(poly3, exps3, len3)`` using Johnson's heap method (see papers by
    Michael Monagan and Roman Pearce). The function realocates its output, hence
    the double indirection, and returns the length of the product. The function
    assumes the exponent vectors take N words. No aliasing is allowed.

.. function:: void fmpz_mpoly_mul_johnson(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times ``poly3`` using the Johnson heap
    based method. See the numerous papers by Michael Monagan and Roman Pearce.

.. function:: void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Does the same operation as ``fmpz_mpoly_mul_johnson`` but with
    multiple threads.

.. function:: slong _fmpz_mpoly_mul_array(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits)

    Set ``(poly1, exp1, alloc)`` to ``(poly2, exps2, len2)`` times
    ``(poly3, exps3, len3)`` by accumulating coefficients in a big, dense
    array. The function realocates its output, hence the double indirection, and
    returns the length of the product. The array ``mults`` is a list of bases
    to be used in encoding the array indices from the exponents. They should
    exceed the maximum exponent for each field of the exponent vectors of the
    output. The output exponent vectors will be packed with fields of the given
    number of bits. The number of variables is given by ``num``. No aliasing
    is allowed.

.. function:: int fmpz_mpoly_mul_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times ``poly3`` using a big array to
    accumulate coefficients. If the array will be larger than some internally
    set parameter, the function fails silently and returns 0 so that some other
    method may be called. This function is most efficient on semi-sparse inputs.

.. function:: int fmpz_mpoly_mul_dense(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, set ``poly1`` to ``poly2`` times
    ``poly3`` using a Kronecker substitution and ``fmpz_poly_mul``.


Powering
----------------------------------------------------------------------


.. function:: slong _fmpz_mpoly_pow_fps(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, slong k, slong N)

    Set ``(poly2, exp1, alloc)`` ``(poly2, exp2, len2)`` raised to the
    power of `k`. The function reallocates its output, hence the double
    indirection. Assumes that exponents vectors each take `N` words. Uses the
    FPS algorithm of Monagan and Pearce. No aliasing is allowed. Assumes
    ``len2 > 1``.

.. function:: void fmpz_mpoly_pow_fps(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong k, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` raised to the `k`-th power, using the
    Monagan and Pearce FPS algorithm. It is assumed that `k \geq 0`.


Divisibility testing
----------------------------------------------------------------------


.. function:: slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits)

    Use dense array exact division to set ``(poly1, exp1, alloc)`` to
    ``(poly2, exp3, len2)`` divided by ``(poly3, exp3, len3)`` in
    ``num`` variables, given a list of multipliers to tightly pack exponents
    and a number of bits for the fields of the exponents of the result. The
    array "mults" is a list of bases to be used in encoding the array indices
    from the exponents. The function reallocates its output, hence the double
    indirection and returns the length of its output if the quotient is exact,
    or zero if not. It is assumed that ``poly2`` is not zero. No aliasing is
    allowed.

.. function:: int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by ``poly3``, using a big dense
    array to accumulate coefficients and return 1 if the quotient is exact.
    Otherwise, return 0 if the quotient is not exact. If the array will be
    larger than some internally set parameter, the function fails silently and
    returns `-1` so that some other method may be called. This function is most
    efficient on dense inputs. Note that the function 
    ``fmpz_mpoly_div_monagan_pearce`` below may be much faster if the
    quotient is known to be exact.

.. function:: slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N)

    Set ``(poly1, exp1, alloc)`` to ``(poly2, exp3, len2)`` divided by
    ``(poly3, exp3, len3)`` and return 1 if the quotient is exact. Otherwise
    return 0. The function assumes exponent vectors that each fit in `N` words,
    and are packed into fields of the given number of bits. Assumes input polys
    are nonzero. Implements "Polynomial division using dynamic arrays, heaps
    and packed exponents" by Michael Monagan and Roman Pearce. No aliasing is
    allowed.

.. function:: int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by ``poly3`` and return 1 if
    the quotient is exact. Otherwise return 0. The function uses the algorithm
    of Michael Monagan and Roman Pearce. Note that the function
    ``fmpz_mpoly_div_monagan_pearce`` below may be much faster if the
    quotient is known to be exact.


Division
----------------------------------------------------------------------


.. function:: slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq, ulong ** expq, slong * allocq, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N)

    Set ``(polyq, expq, allocq)`` to the quotient of
    ``(poly2, exp2, len2)`` by ``(poly3, exp3, len3)`` discarding
    remainder (with notional remainder coefficients reduced modulo the leading
    coefficient of ``(poly3, exp3, len3)``), and return the length of the
    quotient. The function reallocates its output, hence the double
    indirection. The function assumes the exponent vectors all fit in `N`
    words. The exponent vectors are assumed to have fields with the given
    number of bits. Assumes input polynomials are nonzero. Implements
    "Polynomial division using dynamic arrays, heaps and packed exponents" by
    Michael Monagan and Roman Pearce. No aliasing is allowed.

.. function:: void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t polyq, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``polyq`` to the quotient of ``poly2`` by ``poly3``,
    discarding the remainder (with notional remainder coefficients reduced
    modulo the leading coefficient of ``poly3``). Implements "Polynomial
    division using dynamic arrays, heaps and packed exponents" by Michael
    Monagan and Roman Pearce. This function is exceptionally efficient if the
    division is known to be exact.

.. function:: slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N)

    Set ``(polyq, expq, allocq)`` and ``(polyr, expr, allocr)`` to the
    quotient and remainder of ``(poly2, exp2, len2)`` by
    ``(poly3, exp3, len3)`` (with remainder coefficients reduced modulo the
    leading coefficient of ``(poly3, exp3, len3)``), and return the length
    of the quotient. The function reallocates its outputs, hence the double
    indirection. The function assumes the exponent vectors all fit in `N`
    words. The exponent vectors are assumed to have fields with the given
    number of bits. Assumes input polynomials are nonzero. Implements
    "Polynomial division using dynamic arrays, heaps and packed exponents" by
    Michael Monagan and Roman Pearce. No aliasing is allowed.

.. function:: void fmpz_mpoly_divrem_monagan_pearce(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``polyq`` and ``polyr`` to the quotient and remainder of
    ``poly2`` divided by ``poly3``, (with remainder coefficients reduced
    modulo the leading coefficient of ``poly3``). Implements "Polynomial
    division using dynamic arrays, heaps and packed exponents" by Michael
    Monagan and Roman Pearce.

.. function:: slong _fmpz_mpoly_divrem_array(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits)

    Use dense array division to set ``(polyq, expq, allocq)`` and
    ``(polyr, expr, allocr)`` to the quotient and remainder of
    ``(poly2, exp2, len2)`` divided by ``(poly3, exp3, len3)`` in
    ``num`` variables, given a list of multipliers to tightly pack
    exponents and a number of bits for the fields of the exponents of the
    result. The function reallocates its outputs, hence the double indirection.
    The array ``mults`` is a list of bases to be used in encoding the array
    indices from the exponents. The function returns the length of the
    quotient. It is assumed that the input polynomials are not zero. No
    aliasing is allowed.

.. function:: int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``polyq`` and ``polyr`` to the quotient and remainder of
    ``poly2`` divided by ``poly3``, (with remainder coefficients reduced
    modulo the leading coefficient of ``poly3``). The function is
    implemented using dense arrays, and is efficient when the inputs are fairly
    dense. If the array will be larger than some internally set parameter, the
    function silently returns 0 so that another function can be called,
    otherwise it returns 1.

.. function:: void fmpz_mpoly_quasidivrem_heap(fmpz_t scale, fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``scale``, ``q`` and ``r`` so that
    ``scale*poly2 = q*poly3 + r`` and no monomial in ``r`` is divisible
    by the leading monomial of ``poly3``, where ``scale`` is positive
    and as small as possible. This function throws an execption if
    ``poly3`` is zero or if an exponent overflow occurs.


Reduction
----------------------------------------------------------------------


.. function:: slong _fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3, ulong * const * exp3, slong len, slong N, slong bits, const fmpz_mpoly_ctx_t ctx)

    This function is as per ``_fmpz_mpoly_divrem_monagan_pearce`` except
    that it takes an array of divisor polynomials ``poly3`` and an array of
    repacked exponent arrays ``exp3``, which may alias the exponent arrays
    of ``poly3``, and it returns an array of quotient polynomials
    ``polyq``. The number of divisor (and hence quotient) polynomials, is
    given by ``len``. The function computes polynomials `q_i` such that
    `r = a - \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where the `q_i` are the
    quotient polynomials and the `b_i` are the divisor polynomials.

.. function:: void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len, const fmpz_mpoly_ctx_t ctx)

    This function is as per ``fmpz_mpoly_divrem_monagan_pearce`` except
    that it takes an array of divisor polynomials ``poly3``, and it returns
    an array of quotient polynomials ``q``. The number of divisor (and hence
    quotient) polynomials, is given by ``len``. The function computes
    polynomials `q_i = q[i]` such that ``poly2`` is
    `r + \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where `b_i =` ``poly3[i]``.


Differentiation/Integration
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_derivative(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong idx, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to the derivative of ``poly2`` with respect to the
    variable of index ``idx``. This function cannot fail.

.. function:: void fmpz_mpoly_integral(fmpz_mpoly_t poly1, fmpz_t scale, const fmpz_mpoly_t poly2, slong idx, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` and ``scale`` so that ``poly1`` is an integral of
    ``poly2*scale`` with respect to the variable of index ``idx``,
    where ``scale`` is positive and as small as possible. This function
    throws an exception upon exponent overflow.


Evaluation
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A, fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)

    Set ``ev`` to the evaluation of ``poly`` where the variables are
    replaced by the corresponding elements of the array ``vals``. This
    function uses a tree method on the variable of largest degree.

.. function:: void fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of
    index ``var`` is replaced by ``val``.

.. function:: void fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctxB)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.

.. function:: void fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``. Both
    ``A`` and the elements of ``C`` have context object
    ``ctxAC``, while ``B`` has context object ``ctxB``. Neither of
    ``A`` and ``B`` is allowed to alias any other polynomial.


Greatest Common Divisor
----------------------------------------------------------------------


.. function:: void fmpz_mpoly_term_content(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Sets ``poly1`` to the GCD of the terms of ``poly2``.


.. function:: int fmpz_mpoly_gcd_prs(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, used psuedo-remainder sequences to set 
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` has positive leading term.

.. function:: int fmpz_mpoly_gcd_brown(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, used Brown's dense modular algorithm to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` has positive leading term.

.. function:: int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, used a modular algorithm with Zippel's sparse
    interpolation to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` has positive leading term.

.. function:: int fmpz_mpoly_resultant(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, slong var, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, set ``poly1`` to the resultant of 
    ``poly2`` and ``poly3`` with respect to the variable of
    index ``var``.

.. function:: int fmpz_mpoly_discriminant(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, set ``poly1`` to the discriminant of
    ``poly2`` with respect to the variable of index ``var``.


Univariates
----------------------------------------------------------------------



.. function:: void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)

    Initialize ``poly``.

.. function:: void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t poly, const fmpz_mpoly_ctx_t ctx)

    Free all memory used by ``poly``.

.. function:: void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t poly1, fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)

    Swap ``poly1`` and ``poly2``.

.. function:: void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t poly, slong length, const fmpz_mpoly_ctx_t ctx)

    Make sure that ``poly`` has space for at least ``length`` terms.

.. function:: int fmpz_mpoly_to_univar(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)

    If return is nonzero, broke up ``poly2`` as a polynomial
    in the variable of index ``var``
    with multivariate coefficients in the other variables, and stored the result
    in ``poly1``. The return is zero if and only if the degree of
    ``poly2`` with respect to the variable of index ``var`` is greater
    or equal to ``2^(FLINT_BITS-1)``.

.. function:: void fmpz_mpoly_from_univar(fmpz_mpoly_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)

    Reverse the operation performed by ``fmpz_mpoly_to_univar``. This
    function is currently undefined if the coefficients of ``poly2``
    themselves depend on the main variable in ``poly2``. 

.. function:: int fmpz_mpoly_univar_equal(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)

    Return 1 if ``poly1`` and ``poly2`` are equal, otherwise return 0.

.. function:: void fmpz_mpoly_univar_add(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus ``poly3``.

.. function:: int fmpz_mpoly_univar_mul(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3, const fmpz_mpoly_ctx_t ctx)

    If return is nonzero, set ``poly1`` to ``poly2`` times ``poly3``.

.. function:: void fmpz_mpoly_univar_derivative(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to the derivative of ``poly2`` with respect to
    its main variable.

.. function:: void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * shift1, const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` and ``shift1`` so that `p_1*x^{s_1} = p_2`. The
    shift is included because the ``fmpz_poly_t`` type is a dense type and
    ``fmpz_mpoly_t`` is not. A call to
    ``fmpz_poly_shift_left(poly1, poly1, shift1)``
    will result in ``poly1`` being equal to ``poly2``. This function
    is defined only if ``poly2`` depends solely on the variable
    of index ``var``.

.. function:: void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1, const fmpz_poly_t poly2, slong shift2, slong var, const fmpz_mpoly_ctx_t ctx)

    Reverse the operation performed by ``fmpz_mpoly_to_fmpz_poly``.

.. function:: void _fmpz_mpoly_univar_prem(fmpz_mpoly_univar_t polyA, const fmpz_mpoly_univar_t polyB, fmpz_mpoly_univar_t polyC, const fmpz_mpoly_ctx_t ctx)

    Set ``polyA`` to the pseudo remainder of ``polyA`` and -``polyB``.
    The division is performed with respect to the variable store in
    ``polyB``. An extra polynomial ``polyC`` is needed for workspace.

.. function:: void _fmpz_mpoly_univar_pgcd(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ, const fmpz_mpoly_ctx_t ctx)

    Set ``poly1`` to the last (nonzero) subresultant polynomial of
    ``polyQ`` and ``polyQ``. It is assumed that `\operatorname{deg}(P)
    \ge \operatorname{deg}(Q) \ge 1`.

.. function:: void _fmpz_mpoly_univar_pgcd_ducos(fmpz_mpoly_univar_t poly1, const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ, const fmpz_mpoly_ctx_t ctx)

    Perform the same operation as ``_fmpz_mpoly_univar_pgcd`` using the
    algorithm of Ducos.


Input/Output
----------------------------------------------------------------------


.. function:: char * _fmpz_mpoly_get_str_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, slong n, int deg, int rev, slong N)

    Returns a string (which the user is responsible for cleaning up),
    representing ``(poly, exps, len)`` in `n` variables, exponent fields
    of the given number of bits and exponent vectors taking `N` words each,
    given an array of `n` variable strings, starting with the variable of
    most significance with respect to the ordering. The ordering is
    specified by the values ``deg``, which is set to 1 if the polynomial
    is deglex or degrevlex, and ``rev``, which is set to 1 if the
    polynomial is revlex or degrevlex.

.. function:: char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Return a string (which the user is responsible for cleaning up),
    representing ``poly``, given an array of variable strings, starting
    with the variable of most significance with respect to the ordering. 

.. function:: int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, slong n, int deg, int rev, slong N)

    Print to the given stream, a string representing ``(poly, exps, len)``
    in `n` variables, exponent fields of the given number of bits and exponent
    vectors taking `N` words each, given an array of `n` variable strings,
    starting with the variable of most significance with respect to the
    ordering. The ordering is specified by the values ``deg``, which is set
    to 1 if the polynomial is deglex or degrevlex, and ``rev``, which is set
    to 1 if the polynomial is revlex or degrevlex. The number of characters
    written is returned.

.. function:: int fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Print to the given stream, a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int _fmpz_mpoly_print_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, slong n, int deg, int rev, slong N)

    Print to stdout, a string representing ``(poly, exps, len)``
    in `n` variables, exponent fields of the given number of bits and exponent
    vectors taking `N` words each, given an array of `n` variable strings,
    starting with the variable of most significance with respect to the
    ordering. The ordering is specified by the values ``deg``, which is set
    to 1 if the polynomial is deglex or degrevlex, and ``rev``, which is set
    to 1 if the polynomial is revlex or degrevlex. The number of characters
    written is returned.

.. function:: int fmpz_mpoly_print_pretty(const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Print to the given stream, a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t poly, const char * str, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Sets ``poly`` to the polynomial in the null-terminates string ``str``
    given an array ``x`` of variable strings. If parsing ``str`` fails,
    ``poly`` is set to zero, and ``-1`` is returned. Otherwise, ``0``
    is returned. The operations ``+``, ``-``, ``*``, and ``/`` are
    permitted along with integers and the variables in ``x``. The character
    ``^`` must be immediately followed by the (integer) exponent. If any
    division is not exact, parsing fails.


Random generation
----------------------------------------------------------------------

.. function:: void fmpz_mpoly_randtest_bound(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong exp_bound, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to ``n_randint(state, exp_bound)``.

.. function:: void fmpz_mpoly_randtest_bounds(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong * exp_bounds, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index ``i`` are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, mp_limb_t exp_bits, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and exponents whose packed form does not exceed the given bit count.

    The parameter ``coeff_bits`` to the three functions ``fmpz_mpoly_randtest_{bound|bounds|bits}`` is merely a suggestion for the approximate bit count of the resulting signed coefficients.
    The function ``fmpz_mpoly_max_bits`` will give the exact bit count of the result.
