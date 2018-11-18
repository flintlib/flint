.. _fmpq-mpoly:

**fmpq_mpoly.h** -- multivariate polynomials over the rational numbers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpq_mpoly_ctx_struct

.. type:: fmpq_mpoly_ctx_t

    Description.

.. type:: fmpq_mpoly_struct

.. type:: fmpq_mpoly_t

    Description.


Context object
--------------------------------------------------------------------------------


    An array of type ``ulong *`` or ``fmpz **`` is used to communicate
    exponent vectors. These exponent vectors must have length equal to the
    number of variables in the polynmial ring.
    The element of this exponent vector at index `0`
    corresponds to the most significant variable in the monomial ordering.
    For example, if the polynomial is `7*x^2*y+8*y*z+9` and the variables are
    ordered so that `x>y>z`, the degree function will return `{2,1,1}`.
    Similarly, the exponent vector of the `0`-index term in this polynomial is
    `{2,1,0}`, while the `2`-index term has exponent vector `{0,0,0}` and
    coefficient `9`.

.. function:: void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

    Initialise a context object for a polynomial ring with the given number of
    variables and the given ordering. The possibilities for the ordering are
    ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong fmpq_mpoly_ctx_nvars(fmpq_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)

    Release up any space allocated by an ``fmpq_mpoly_ctx_t``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_init(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Initialise ``poly`` for use, given an initialised context object.
    Its value is set to 0.

.. function:: void fmpq_mpoly_init2(fmpq_mpoly_t poly, slong alloc, const fmpq_mpoly_ctx_t ctx)

    Initialise ``poly`` for use, with space for at least
    ``alloc`` terms, given an initialised context. Its value is set to 0.

.. function:: void fmpq_mpoly_realloc(fmpq_mpoly_t poly, slong len, const fmpq_mpoly_ctx_t ctx)

    Reallocate ``poly`` to have space for ``len`` terms. 
    Assumes the current length of the polynomial is not greater than
    ``len``.

.. function:: void fmpq_mpoly_fit_length(fmpq_mpoly_t poly, slong len, const fmpq_mpoly_ctx_t ctx)

    Reallocate ``poly`` to have space for at least
    ``len`` terms. No truncation is performed if ``len`` is less than
    the currently allocated number of terms; the allocated space can only grow.

.. function:: void fmpq_mpoly_fit_bits(fmpq_mpoly_t poly, slong bits, const fmpq_mpoly_ctx_t ctx)

    Reallocate the polynomial to have space for exponent fields of the given
    number of bits. This function can increase the number of bits only.

.. function:: void fmpq_mpoly_clear(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Release any space allocated an ``fmpq_mpoly_t``.


Degrees
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Return 1 if the degrees of poly with respect to each variable fit into
    an ``slong``, otherwise return 0.

.. function:: void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``poly`` with respect to each variable.
    If ``poly`` is zero, all degrees are set to ``-1``.

.. function:: slong fmpq_mpoly_degree_si(const fmpq_mpoly_t poly, slong var, const fmpq_mpoly_ctx_t ctx)

    Return the degree of ``poly`` with respect to the variable of index
    ``var``. If ``poly`` is zero, the return is ``-1``.

.. function:: void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``poly`` with respect to each variable.
    If ``poly`` is zero, all degrees are set to ``-1``.

.. function:: void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t poly, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``deg`` to the degree of ``poly`` with respect to the variable
    of index ``var``. If ``poly`` is zero, set ``deg`` to ``-1``.

.. function:: int fmpq_mpoly_totaldegree_fits_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return 1 if the total degree of ``A`` fits into
    an ``slong``, otherwise return 0.

.. function:: slong fmpq_mpoly_totaldegree_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return the total degree of ``A`` assuming it fits into an slong.
    If ``A`` is zero, the return is ``-1``.

.. function:: void fmpq_mpoly_totaldegree_fmpz(fmpz_t tdeg, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, ``tdeg`` is set to ``-1``.


Coefficients
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_denominator(fmpz_t d, const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Set `d` to the denominator of `poly`, the smallest positive integer `d`
    such that `d`*``poly`` has integer coefficients.

.. function:: void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c, const fmpq_mpoly_t poly, const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``poly2`` is a monomial,
    set `c` to the coefficient of the corresponding monomial in ``poly``.
    This function thows if ``poly2`` is not a monomial.

.. function:: void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t poly, const fmpq_t c, const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``poly2`` is a monomial,
    set the coefficient of the corresponding monomial in ``poly`` to `c`.
    This function thows if ``poly2`` is not a monomial.

.. function:: void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t poly, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)

    Set `c` to the coefficient of the monomial with exponent ``exp``.

.. function:: void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t poly, const fmpq_t c, fmpz * const * exp, fmpq_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent ``exp`` to `c`.

.. function:: void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t poly, ulong const * exp, const fmpq_mpoly_ctx_t ctx)

    Set `c` to the coefficient of the monomial with exponent ``exp``.

.. function:: void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t poly, const fmpq_t c, ulong const * exp, fmpq_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent ``exp`` to `c`.


Internal operations
--------------------------------------------------------------------------------


    An ``fmpq_mpoly_t poly`` is an ``fmpz_mpoly_t zpoly`` together with
    an ``fmpq_t content`` representing any content.
    In order to be in canonical form, either ``zpoly`` and ``content``
    must both be zero, or ``zcode`` must be a polynomial with content 1
    and positive leading term. In either the case the representation
    ``poly = content * zpoly`` holds.

.. function:: void fmpq_mpoly_canonicalise(fmpq_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Factor out content from the internal ``fmpz_mpoly_t``. If this
    ``fmpz_mpoly_t`` was in canonical form, the resulting ``poly``
    will be. This function expects the internal ``fmpz_mpoly_t`` to be
    in canonical form. If this is not the case, the container operations
    of ``fmpq_mpoly_sort`` and ``fmpq_mpoly_combine_like_terms`` may
    be useful. All operations other than the container operations 
    produce canonical output and expect their inputs to be canonical.

.. function:: void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Throw if ``poly`` is not in canonical form.


Container operations
--------------------------------------------------------------------------------


    Some of these functions deal with violations of the internal canonical
    representation. A call to ``fmpq_mpoly_sort`` followed by a call to
    ``fmpq_mpoly_combine_like_terms`` should leave a polynomial in
    canonical form. The ``pushterm`` functions run in constant
    average time if the terms pushed have bounded denominators,
    and a term is appened even if the specified coefficient is zero.


.. function:: slong fmpq_mpoly_length(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Return the number of terms stored in ``poly``. If the polynomial
    is in canonical form, this will be the number of nonzero coefficients.


.. function:: void fmpq_mpoly_get_termcoeff_fmpq(fmpq_t x, const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Set `x` to coefficient of index `i`, starting
    with `i = 0` for the term with most significance.

.. function:: void fmpq_mpoly_set_termcoeff_fmpq(fmpq_mpoly_t poly, slong i, const fmpq_t x, const fmpq_mpoly_ctx_t ctx)

    Set the coefficient of index `i` to `x`, starting
    with `i = 0` for the term with most significance.


.. function:: int fmpq_mpoly_termexp_fits_si(const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Return 1 if all entries of the exponent vector of the term of index `i`
    fit into an ``slong``. Otherwise, return 0.

.. function:: int fmpq_mpoly_termexp_fits_ui(const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Return 1 if all entries of the exponent vector of the term of index `i`
    fit into a ``ulong``. Otherwise, return 0.

.. function:: void fmpq_mpoly_get_termexp_ui(ulong * exps, const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Get the exponent vector of the given polynomial with index `i`.

.. function:: void fmpq_mpoly_get_termexp_fmpz(fmpz ** exps, const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Get the exponent vector of the given polynomial with index `i`.

.. function:: void fmpq_mpoly_set_termexp_ui(fmpq_mpoly_t poly, slong i, const ulong * exps, const fmpq_mpoly_ctx_t ctx)

    Set the exponent vector of the given polynomial with index `i`.

.. function:: void fmpq_mpoly_set_termexp_fmpz(fmpq_mpoly_t poly, slong i, fmpz * const * exps, const const fmpq_mpoly_ctx_t ctx)

    Set the exponent vector of the given polynomial with index `i`.


.. function:: void fmpq_mpoly_pushterm_fmpq_fmpz(fmpz_mpoly_t poly, const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_fmpz_fmpz(fmpz_mpoly_t poly, const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_ui_fmpz(fmpz_mpoly_t poly, ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_si_fmpz(fmpz_mpoly_t poly, slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_fmpq_ui(fmpz_mpoly_t poly, const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_fmpz_ui(fmpz_mpoly_t poly, const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_ui_ui(fmpz_mpoly_t poly, ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.

.. function:: void fmpq_mpoly_pushterm_si_ui(fmpz_mpoly_t poly, slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``poly`` with the given coefficient and exponents.


.. function:: void fmpq_mpoly_sort_terms(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Sort the internal ``fmpz_mpoly_t`` into the canonical ordering
    dictated by the ordering in ``ctx``. This function does not combine
    like terms, nor does it delete terms with coefficient zero.
    Even if all terms have distinct exponents, the result may still not be
    in canonical form, because content may not be factored out.

.. function:: void fmpq_mpoly_combine_like_terms(fmpq_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)

    Combine adjacent like terms in the internal ``fmpz_poly`` and then
    factor out content via ``fmpq_mpoly_canonicalise``. If the terms of
    ``poly`` were sorted to begin with, the result will be in canonical form.


Set/Negate
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_set(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
    
    Set ``poly1`` to ``poly2``.

.. function:: void fmpq_mpoly_swap(fmpq_mpoly_t poly1, fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)

    Efficiently swap the contents of the two given polynomials. No copying is
    performed; the swap is accomplished by swapping pointers.

.. function:: void fmpq_mpoly_gen(fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    Set ``poly`` to the `i`-th generator (variable),
    where `i = 0` corresponds to the variable with the most significance
    with respect to the ordering. 

.. function:: void fmpq_mpoly_neg(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
    
    Set ``poly1`` to `-```poly2``.



Constants
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is a constant, else return 0.

.. function:: void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``poly`` is a constant, set `c` to this constant.
    This function throws if ``poly`` is not a constant.

.. function:: void fmpq_mpoly_set_fmpq(fmpq_mpoly_t poly, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly`` to the constant `c`.

.. function:: void fmpq_mpoly_set_fmpz(fmpq_mpoly_t poly, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly`` to the constant `c`.

.. function:: void fmpq_mpoly_set_ui(fmpq_mpoly_t poly, ulong c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly`` to the constant `c`.

.. function:: void fmpq_mpoly_set_si(fmpq_mpoly_t poly, slong c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly`` to the constant `c`.

.. function:: void fmpq_mpoly_zero(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Set ``poly`` to the constant 0.

.. function:: void fmpq_mpoly_one(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Set ``poly`` to the constant 1.


Comparison
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_equal(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly1`` is equal to ``poly2``, else return 0.

.. function:: int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t poly, fmpq_t c, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t poly, fmpz_t c, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpq_mpoly_equal_ui(const fmpq_mpoly_t poly, ulong  c, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpq_mpoly_equal_si(const fmpq_mpoly_t poly, slong  c, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant `c`, else return 0.

.. function:: int fmpq_mpoly_is_zero(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant 0, else return 0.

.. function:: int fmpq_mpoly_is_one(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

    Return 1 if ``poly`` is equal to the constant 1, else return 0.


.. function:: int fmpq_mpoly_is_gen(const fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)

    If `i \ge 0`, return 1 if ``poly`` is equal to the `i`-th generator,
    otherwise return 0. If `i < 0`, return 1 if the polynomial is
    equal to any generator, otherwise return 0.


Basic arithmetic
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_add_fmpq(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpq_mpoly_add_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpq_mpoly_add_ui(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpq_mpoly_add_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` plus `c`.

.. function:: void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpq_mpoly_sub_ui(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpq_mpoly_sub_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

    Set ``poly1`` to ``poly2`` minus `c`.

.. function:: void fmpq_mpoly_add(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus ``poly3``.

.. function:: void fmpq_mpoly_sub(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus ``poly3``.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times `c`.

.. function:: void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times `c`.

.. function:: void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times `c`.

.. function:: void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times `c`.

.. function:: void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by `c`.

.. function:: void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by `c`.

.. function:: void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by `c`.

.. function:: void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by `c`.

.. function:: void fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t poly1, const fmpq_mpoly_ctx_t ctx)

    Divide ``poly1`` by its leading coefficient. An expection is raised if
    ``poly1`` is zero.


Multiplication
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_mul(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times ``poly3``.

.. function:: void fmpq_mpoly_mul_threaded(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times ``poly3`` using multiple threads.


Powering
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_pow_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpz_t k, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` raised to the `k`-th power.
    An expection is raised if `k < 0` or if `k` is large and the polynomial is
    not a monomial with coefficient `\pm1`.

.. function:: void fmpq_mpoly_pow_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong k, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` raised to the `k`-th power.
    An expection is raised if `k < 0`.


Divisibility testing
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_divides(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by ``poly3`` and return 1 if
    the quotient is exact. Otherwise return 0.


Division
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_div(fmpq_mpoly_t polyq, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``polyq`` to the quotient of ``poly2`` by ``poly3``,
    discarding the remainder. An expection is currently raised if ``poly2``
    or ``poly3`` have bit counts greater than ``FLINT_BITS``.

.. function:: void fmpq_mpoly_divrem(fmpq_mpoly_t q, fmpq_mpoly_t r, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``polyq`` and ``polyr`` to the quotient and remainder of
    ``poly2`` divided by ``poly3``. An expection is 
    currently raised if ``poly2``
    or ``poly3`` have bit counts greater than ``FLINT_BITS``.


Reduction
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** q, fmpq_mpoly_t r, const fmpq_mpoly_t poly2, fmpq_mpoly_struct * const * poly3, slong len, const fmpq_mpoly_ctx_t ctx)

    This function is as per ``fmpq_mpoly_divrem`` except
    that it takes an array of divisor polynomials ``poly3``, and it returns
    an array of quotient polynomials ``q``. The number of divisor (and hence
    quotient) polynomials, is given by ``len``. The function computes
    polynomials `q_i = q[i]` such that ``poly2`` is
    `r + \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where `b_i =` ``poly3[i]``.
    An expection is currently raised 
    if any input polynomials have bit counts greater than ``FLINT_BITS``.


Differentiation/Integration
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_derivative(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to the derivative of ``poly2`` with respect to the
    variable of index ``var``.

.. function:: void fmpq_mpoly_integral(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to the integral with the fewest number of terms
    of ``poly2`` with respect to the variable of index ``var``.


Evaluation
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A, fmpq * const * vals, const fmpq_mpoly_ctx_t ctx)

    Set ``ev`` the evaluation of ``A`` where the variables are
    replaced by the corresponding elements of the array ``vals``.

.. function:: void fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_t val, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of
    index ``var`` is replaced by ``val``.

.. function:: void fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A, const fmpq_mpoly_t B, fmpq_poly_struct * const * C, const fmpq_mpoly_ctx_t ctxB)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.

.. function:: void fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A, const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C, const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``. Both
    ``A`` and the elements of ``C`` have context object
    ``ctxAC``, while ``B`` has context object ``ctxB``. Neither of
    ``A`` and ``B`` is allowed to alias any other polynomial.


Greatest Common Divisor
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_gcd(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)

    Set ``poly1`` to the monic GCD of ``poly2`` and ``poly3``, assuming
    the return value is 1. If the return value is 0, the GCD was
    unable to be computed.


Input/Output
--------------------------------------------------------------------------------


.. function:: char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t poly, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Return a string (which the user is responsible for cleaning up),
    representing ``poly``, given an array of variable strings, starting
    with the variable of most significance with respect to the ordering. 

.. function:: int fmpq_mpoly_fprint_pretty(FILE * file, const fmpq_mpoly_t poly, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Print to the given stream a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int fmpq_mpoly_print_pretty(const fmpq_mpoly_t poly, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Print to ``stdout`` a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t poly, const char * str, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Sets ``poly`` to the polynomial in the null-terminated string ``str``
    given an array ``x`` of variable strings. If parsing ``str`` fails,
    ``poly`` is set to zero, and ``-1`` is returned. Otherwise, ``0``
    is returned. The operations ``+``, ``-``, ``*``, and ``/`` are
    permitted along with integers and the variables in ``x``. The character
    ``^`` must be immediately followed by the (integer) exponent. If any
    division is not exact, parsing fails.


Random generation
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_randtest_bound(fmpq_mpoly_t poly, flint_rand_t state, slong length, mp_limb_t coeff_bits, slong exp_bound, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with
    length up to the given length,
    exponents in the range ``[0, exp_bound - 1]``, and with
    rational coefficients of the given number of bits.
    The exponents of each variable are generated by calls to
    ``n_randint(state, exp_bound)``.

.. function:: void fmpq_mpoly_randtest_bound(fmpq_mpoly_t poly, flint_rand_t state, slong length, mp_limb_t coeff_bits, slong exp_bound, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with
    length up to the given length,
    exponents in the range ``[0, exp_bounds[i] - 1]``, and with
    rational coefficients of the given number of bits.
    The exponents of the variable of index `i` are generated by calls to
    ``n_randint(state, exp_bounds[i])``.

.. function:: void fmpq_mpoly_randtest_bits(fmpq_mpoly_t poly, flint_rand_t state, slong length, mp_limb_t coeff_bits, mp_limb_t exp_bits, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length, exponents
    whose packed form does not exceed the given bit count, and with rational
    coefficients of the given number of bits.
