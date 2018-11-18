.. _nmod-mpoly:

**nmod_mpoly.h** -- multivariate polynomials over integers mod n (word-size n)
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: nmod_mpoly_ctx_struct

.. type:: nmod_mpoly_ctx_t

    Description.

.. type:: nmod_mpoly_struct

.. type:: nmod_mpoly_t

    Description.


Context object
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, mp_limb_t n)

    Initialise a context object for a polynomial ring with the given number of
    variables and the given ordering.  It will have coefficients modulo~`n`.
    The possibilities for the ordering are
    ``ORD_LEX``, ``ORD_REVLEX``, ``ORD_DEGLEX`` and
    ``ORD_DEGREVLEX``.

.. function:: void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)

    Release any space allocated by an ``nmod_mpoly_ctx_t``.

.. function:: void nmodf_ctx_init(nmodf_ctx_t ctx, mp_limb_t n)

    Initialise a context object for doing arithmetic modulo~`n`.

.. function:: void nmodf_ctx_clear(nmodf_ctx_t ctx)

    Release any space allocated by an ``nmodf_ctx_t``.


Memory management
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_init(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Initialise an ``nmod_mpoly_t`` for use, given an initialised context
    object.

.. function:: void nmod_mpoly_init2(nmod_mpoly_t poly, slong alloc, const nmod_mpoly_ctx_t ctx)

    Initialise an ``nmod_mpoly_t`` for use, with space for at least
    ``alloc`` terms, given an initialised context. By default, fields of 8
    bits are allocated for the exponents in each exponent vector.

.. function:: void nmod_mpoly_realloc(nmod_mpoly_t poly, slong len, const nmod_mpoly_ctx_t ctx)

    Reallocate an ``nmod_mpoly_t`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than
    ``len``.

.. function:: void nmod_mpoly_clear(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Release any space allocated for an ``nmod_mpoly_t``.

.. function:: void nmod_mpoly_truncate(nmod_mpoly_t poly, slong newlen, const nmod_mpoly_ctx_t ctx)

    If the given polynomial is larger than the given number of terms, truncate
    to that number of terms.

.. function:: void nmod_mpoly_fit_bits(nmod_mpoly_t poly, slong bits, const nmod_mpoly_ctx_t ctx)

    Reallocate the polynomial to have space for exponent fields of the given
    number of bits. The number of bits must be at least 8 and at most
    FLINT_BITS. This function can increase the number of bits only.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_gen(nmod_mpoly_t poly, slong i, const nmod_mpoly_ctx_t ctx)

    Set the given ``nmod_mpoly_t`` to the `i`-th generator (variable),
    where `i = 0` corresponds to the variable with the most significance
    with respect to the ordering. 

.. function:: void nmod_mpoly_set_ui(nmod_mpoly_t poly, ulong c, const nmod_mpoly_ctx_t ctx)

    Set the given ``nmod_mpoly_t`` to the constant polynomial corresponding
    to the unsigned integer `c`.

.. function:: int nmod_mpoly_equal_ui(const nmod_mpoly_t poly, ulong c, const nmod_mpoly_ctx_t ctx)

    Return 1 if the given ``nmod_mpoly_t`` is equal to the constant
    polynomial represented by the unsigned integer `c`, else return 0.

.. function:: void nmod_mpoly_swap(nmod_mpoly_t poly1, nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)

    Efficiently swap the contents of the two given polynomials. No copying is
    performed; the swap is accomplished by swapping pointers.

.. function:: void nmod_mpoly_zero(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Set the given ``nmod_mpoly_t`` to the zero polynomial.

.. function:: void nmod_mpoly_one(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Set the given ``nmod_mpoly_t`` to the constant polynomial with value 1.

.. function:: int nmod_mpoly_is_zero(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Return 1 if the given ``nmod_mpoly_t`` is equal to the zero polynomial,
    else return 0.

.. function:: int nmod_mpoly_is_one(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Return 1 if the given ``nmod_mpoly_t`` is equal to the constant
    polynomial with coefficient 1, else return 0.

.. function:: int nmod_mpoly_is_gen(const nmod_mpoly_t poly, slong i, const nmod_mpoly_ctx_t ctx)

    Return 1 if the given ``nmod_mpoly_t`` is equal to the `i`-th generator
    of the polynomial ring, otherwise return 0. The generator corresponding to
    the variable with the most significance corresponds to `i = 0`. If
    `i < 0` the function returns 1 if the polynomial is equal to any generator
    of the polynomial ring, otherwise it returns 0.

.. function:: void nmod_mpoly_set_term_ui_ui(nmod_mpoly_t poly, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)

    Set the term of ``poly`` with the given monomial to the given
    coefficient. The monomial is specified as a vector of exponents with as
    many variables as the polynomial. The most significant variable with
    respect to the ordering is at index 0 of the vector. If a term with that
    monomial already exists in the polynomial, it is overwritten. The term is
    removed if the coefficient is zero. If a term with that monomial doesn't
    exist, one is inserted at the appropriate position.

.. function:: void nmod_mpoly_set_term_ui_fmpz(nmod_mpoly_t poly, ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)

    Set the term of ``poly`` with the given monomial to the given
    coefficient. The monomial is specified as a vector of exponents with as
    many variables as the polynomial. The most significant variable with
    respect to the ordering is at index 0 of the vector. If a term with that
    monomial already exists in the polynomial, it is overwritten. The term is
    removed if the coefficient is zero. If a term with that monomial doesn't
    exist, one is inserted at the appropriate position.

.. function:: ulong nmod_mpoly_get_term_ui_ui(const fmpz_mpoly_t poly, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    Get the coefficient of the term of ``poly`` with the given monomial.
    The monomial is specified as a vector of exponents with as
    many variables as the polynomial. The most significant variable with
    respect to the ordering is at index 0 of the vector. If no term with that
    monomial exists in the polynomial, zero is returned.

.. function:: ulong nmod_mpoly_get_term_ui_fmpz(const fmpz_mpoly_t poly, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    Get the coefficient of the term of ``poly`` with the given monomial.
    The monomial is specified as a vector of exponents with as
    many variables as the polynomial. The most significant variable with
    respect to the ordering is at index 0 of the vector. If no term with that
    monomial exists in the polynomial, zero is returned.


Degrees
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Return 1 if the degrees of poly with respect to each variable fit into
    an ``slong``, otherwise return 0.

.. function:: void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``poly`` with respect to each variable.
    If ``poly`` is zero, all degrees are set to ``-1``.

.. function:: slong nmod_mpoly_degree_si(const nmod_mpoly_t poly, slong var, const nmod_mpoly_ctx_t ctx)

    Return the degree of ``poly`` with respect to the variable of index
    ``var``. If ``poly`` is zero, the return is ``-1``.

.. function:: void nmod_mpoly_degrees_fmpz(fmpz ** degs, const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``poly`` with respect to each variable.
    If ``poly`` is zero, all degrees are set to ``-1``.

.. function:: void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t poly, slong var, const nmod_mpoly_ctx_t ctx)

    Set ``deg`` to the degree of ``poly`` with respect to the variable
    of index ``var``. If ``poly`` is zero, set ``deg`` to ``-1``.

.. function:: int nmod_mpoly_totaldegree_fits_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return 1 if the total degree of ``A`` fits into
    an ``slong``, otherwise return 0.

.. function:: slong nmod_mpoly_totaldegree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return the total degree of ``A`` assuming it fits into an slong.
    If ``A`` is zero, the return is ``-1``.

.. function:: void nmod_mpoly_totaldegree_fmpz(fmpz_t tdeg, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, ``tdeg`` is set to ``-1``.


Set and negate
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_set(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)
    
    Set ``poly1`` to ``poly2``.

.. function:: void nmod_mpoly_neg(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)
    
    Set ``poly1`` to `-```poly2``.


Comparison
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_equal(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)

    Return 1 if ``poly1`` is equal to ``poly2``, else return 0.


Basic arithmetic
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_add_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus the constant polynomial given by the
    unsigned integer `c`.

.. function:: void nmod_mpoly_sub_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus the constant polynomial given by the
    unsigned integer `c`.

.. function:: void nmod_mpoly_add(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` plus ``poly3``.

.. function:: void nmod_mpoly_sub(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` minus ``poly3``.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times the unsigned integer `c`.


Multiplication
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_mul_johnson(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` times ``poly3`` using the Johnson heap
    based method. See the numerous papers by Michael Monagan and Roman Pearce.
    This function throws an exception upon exponent overflow.


.. function:: void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Does the same operation as ``nmod_mpoly_mul_johnson`` but with
    multiple threads. This function throws an exception upon exponent overflow.


Powering
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_pow(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, ulong k, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` raised to the `k`-th power.

.. function:: void nmod_mpoly_pow_fmpz(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, ulong k, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` raised to the `k`-th power. This function
    will fail if k is bigger than a ``ulong`` and ``poly1`` has more than
    one term.



Divisibility testing
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to ``poly2`` divided by ``poly3`` and return 1 if
    the quotient is exact. Otherwise return 0. The function uses the algorithm
    of Michael Monagan and Roman Pearce. Note that the function
    ``nmod_mpoly_div_monagan_pearce`` below may be much faster if the
    quotient is known to be exact.


Division
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t polyq, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``polyq`` to the quotient of ``poly2`` by ``poly3``,
    discarding the remainder (with notional remainder coefficients reduced
    modulo the leading coefficient of ``poly3``). Implements "Polynomial
    division using dynamic arrays, heaps and packed exponents" by Michael
    Monagan and Roman Pearce. This function is exceptionally efficient if the
    division is known to be exact.

.. function:: void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    Set ``polyq`` and ``polyr`` to the quotient and remainder of
    ``poly2`` divided by ``poly3``, (with remainder coefficients reduced
    modulo the leading coefficient of ``poly3``). Implements "Polynomial
    division using dynamic arrays, heaps and packed exponents" by Michael
    Monagan and Roman Pearce.


Reduction
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** q, nmod_mpoly_t r, const nmod_mpoly_t poly2, nmod_mpoly_struct * const * poly3, slong len, const nmod_mpoly_ctx_t ctx)

    This function is as per ``nmod_mpoly_divrem_monagan_pearce`` except
    that it takes an array of divisor polynomials ``poly3``, and it returns
    an array of quotient polynomials ``q``. The number of divisor (and hence
    quotient) polynomials, is given by ``len``. The function computes
    polynomials `q_i = q[i]` such that ``poly2`` is
    `r + \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where `b_i =` ``poly3[i]``.


Differentiation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_derivative(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, slong idx, const nmod_mpoly_ctx_t ctx)

    Set ``poly1`` to the derivative of ``poly2`` with respect to the
    variable of index ``idx``. This function cannot fail.


Evaluation
--------------------------------------------------------------------------------


.. function:: ulong nmod_mpoly_evaluate_all_ui(nmod_mpoly_t A, const ulong * vals, const nmod_mpoly_ctx_t ctx)

    Return the evaluation of ``A`` where the variables are
    replaced by the corresponding elements of the array ``vals``.

.. function:: void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, ulong val, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of
    index ``var`` is replaced by ``val``.

.. function:: void nmod_mpoly_compose_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``. Both
    ``A`` and the elements of ``C`` should have modulus matching that of
    ``ctxB``, which is the context object of ``B``.

.. function:: void nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are
    replaced by the corresponding elements of the array ``C``. Both
    ``A`` and the elements of ``C`` have context object
    ``ctxAC``, while ``B`` has context object ``ctxB``. Neither of
    ``A`` and ``B`` is allowed to alias any other polynomial.


Greatest Common Divisor
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_gcd_brown(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    If the return is nonzero, used Brown's dense modular algorithm to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` is monic.

.. function:: int nmod_mpoly_gcd_zippel(nmod_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, used a modular algorithm with Zippel's sparse
    interpolation to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` is monic.


Input/Output
--------------------------------------------------------------------------------


.. function:: char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)

    Return a string (which the user is responsible for cleaning up),
    representing ``poly``, given an array of variable strings, starting
    with the variable of most significance with respect to the ordering. 

.. function:: int _nmod_mpoly_fprint_pretty(FILE * file, const nmod * poly, const ulong * exps, slong len, const char ** x, slong bits, slong n, int deg, int rev, slong N)

    Print to the given stream, a string representing ``(poly, exps, len)``
    in `n` variables, exponent fields of the given number of bits and exponent
    vectors taking `N` words each, given an array of `n` variable strings,
    starting with the variable of most significance with respect to the
    ordering. The ordering is specified by the values ``deg``, which is set
    to 1 if the polynomial is deglex or degrevlex, and ``rev``, which is set
    to 1 if the polynomial is revlex or degrevlex. The number of characters
    written is returned.

.. function:: int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)

    Print to the given stream, a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int nmod_mpoly_print_pretty(const nmod_mpoly_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)

    Print to the given stream, a string representing ``poly``, given an
    array of variable strings, starting with the variable of most
    significance with respect to the ordering. The number of characters
    written is returned.

.. function:: int nmod_mpoly_set_str_pretty(nmod_mpoly_t poly, const char * str, const char ** x, const nmod_mpoly_ctx_t ctx)

    Sets ``poly`` to the polynomial in the null-terminates string ``str``
    given an array ``x`` of variable strings. If parsing ``str`` fails,
    ``poly`` is set to zero, and ``-1`` is returned. Otherwise, ``0``
    is returned. The operations ``+``, ``-``, ``*``, and ``/`` are
    permitted along with integers and the variables in ``x``. The character
    ``^`` must be immediately followed by the (integer) exponent. If any
    division is not exact, parsing fails.


Random generation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_randtest_bound(nmod_mpoly_t poly, flint_rand_t state, slong length, slong exp_bound, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with
    length up to the given length and
    exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to
    ``n_randint(state, exp_bound)``.

.. function:: void nmod_mpoly_randtest_bound(nmod_mpoly_t poly, flint_rand_t state, slong length, slong exp_bound, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with
    length up to the given length and
    exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index `i` are generated by calls to
    ``n_randint(state, exp_bounds[i])``.

.. function:: void nmod_mpoly_randtest_bits(nmod_mpoly_t poly, flint_rand_t state, slong length, mp_limb_t exp_bits, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and
    exponents whose packed form does not exceed the given bit count.
