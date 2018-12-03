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

    Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    It will have coefficients modulo `n`. Setting `n = 0` or `n = 1` will currently give undefined behavior.
    The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong nmod_mpoly_ctx_nvars(nmod_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: ordering_t nmod_mpoly_ctx_ord(const nmod_mpoly_ctx_t ctx)

    Return the ordering used to initialize the context.

.. function:: mp_limb_t nmod_mpoly_ctx_modulus(const nmod_mpoly_ctx_t ctx)

    Return the modulus used to initialize the context.

.. function:: void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)

    Release any space allocated by an ``nmod_mpoly_ctx_t``.



Memory management
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    By default 8 bits are allocated for the exponent widths.

.. function:: void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms, and 8 bits are allocated for the exponents.

.. function:: void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms, and ``bits`` bits are allocated for the exponents.

.. function:: void nmod_mpoly_fit_length(nmod_mpoly_t A, slong len, const nmod_mpoly_ctx_t ctx)

    Ensure that ``A`` has space for at least ``len`` terms.

.. function:: void nmod_mpoly_fit_bits(nmod_mpoly_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)

    Ensure that the exponent fields of ``A`` have at least ``bits`` bits.

.. function:: void nmod_mpoly_realloc(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)

    Reallocate ``A`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than ``alloc``.

.. function:: void nmod_mpoly_clear(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Release any space allocated for ``A``.


Input/Output
--------------------------------------------------------------------------------

    The variable strings in ``x`` start with the variable of most significance at index ``0``. If ``x`` is ``NULL``, the variables are named ``x1``, ``x2``, ect.

.. function:: char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing ``A``, given an array of variable strings ``x``.

.. function:: int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``file``.

.. function:: int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``stdout``.

.. function:: int nmod_mpoly_set_str_pretty(nmod_mpoly_t A, const char * str, const char ** x, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the polynomial in the null-terminates string ``str`` given an array ``x`` of variable strings.
    If parsing ``str`` fails, ``A`` is set to zero, and ``-1`` is returned. Otherwise, ``0``  is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in ``x``. The character ``^`` must be immediately followed by the (integer) exponent.
    If any division is not exact, parsing fails.


Basic manipulation
--------------------------------------------------------------------------------

.. function:: void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the variable of index ``var``, where ``var = 0`` corresponds to the variable with the most significance with respect to the ordering. 

.. function:: int nmod_mpoly_is_gen(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

    If `var \ge 0`, return ``1`` if ``A`` is equal to the `var`-th generator, otherwise return ``0``.
    If `var < 0`, return ``1`` if the polynomial is equal to any generator, otherwise return ``0``.

.. function:: void nmod_mpoly_set(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    
    Set ``A`` to ``B``.

.. function:: int nmod_mpoly_equal(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to ``B``, else return ``0``.

.. function:: void nmod_mpoly_swap(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Efficiently swap ``A`` and ``B``.


Constants
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_is_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is a constant, else return ``0``.

.. function:: ulong nmod_mpoly_get_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Assuming that ``A`` is a constant, return this constant.
    This function throws if ``A`` is not a constant.

.. function:: void nmod_mpoly_set_ui(nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``c``.

.. function:: void nmod_mpoly_zero(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``0``.

.. function:: void fmpz_mpoly_one(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``1``.

.. function:: int nmod_mpoly_equal_ui(const nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``c``, else return ``0``.

.. function:: int nmod_mpoly_is_zero(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``0``, else return ``0``.

.. function:: int nmod_mpoly_is_one(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``1``, else return ``0``.


Degrees
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if the degrees of ``A`` with respect to each variable fit into an ``slong``, otherwise return ``0``.

.. function:: void nmod_mpoly_degrees_fmpz(fmpz ** degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

.. function:: void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``A`` with respect to each variable.
    If ``A`` is zero, all degrees are set to ``-1``.

.. function:: void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

.. function:: slong nmod_mpoly_degree_si(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

    Either return or set ``deg`` to the degree of ``A`` with respect to the variable of index ``var``.
    If ``A`` is zero, the degree is defined to be ``-1``.

.. function:: int nmod_mpoly_total_degree_fits_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if the total degree of ``A`` fits into an ``slong``, otherwise return ``0``.

.. function:: void nmod_mpoly_total_degree_fmpz(fmpz_t tdeg, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

.. function:: slong nmod_mpoly_total_degree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Either return or set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, the total degree is defined to be ``-1``.


Coefficients
--------------------------------------------------------------------------------


.. function:: ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, return the coefficient of the corresponding monomial in ``A``.
    This function thows if ``M`` is not a monomial.

.. function:: void nmod_mpoly_set_coeff_ui_monomial(nmod_mpoly_t A, ulong c, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set the coefficient of the corresponding monomial in ``A`` to ``c``.
    This function thows if ``M`` is not a monomial.

.. function:: ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

.. function:: ulong nmod_mpoly_get_coeff_ui_ui(const nmod_mpoly_t A, ulong const * exp, const nmod_mpoly_ctx_t ctx)

    Return the coefficient of the monomial with exponent ``exp``.

.. function:: void nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, nmod_mpoly_ctx_t ctx)

.. function:: void nmod_mpoly_set_coeff_ui_ui(nmod_mpoly_t A, ulong c, ulong const * exp, nmod_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent ``exp`` to `c`.


Container operations
--------------------------------------------------------------------------------

    These functions deal with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.

.. function:: int nmod_mpoly_is_canonical(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is in canonical form. Otherwise, return ``0``.
    To be in canonical form, all of the terms must have nonzero coefficients, and the terms must be sorted from greatest to least.

.. function:: slong nmod_mpoly_length(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Return the number of terms in ``A``.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void nmod_mpoly_resize(nmod_mpoly_t A, slong new_length, const nmod_mpoly_ctx_t ctx)

    Set the length of ``A`` to ``new_length``.
    Terms are either deleted from the end, or new zero terms are appended.

.. function:: ulong nmod_mpoly_get_term_coeff_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

    Return the coefficient of the term of index ``i``.

.. function:: void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c, const nmod_mpoly_ctx_t ctx)

    Set the coefficient of the term of index ``i`` to ``c``.

.. function:: int nmod_mpoly_term_exp_fits_si(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

.. function:: int nmod_mpoly_term_exp_fits_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

    Return ``1`` if all entries of the exponent vector of the term of index `i` fit into an ``slong`` (resp. a ``ulong). Otherwise, return ``0``.

.. function:: void nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

.. function:: void nmod_mpoly_get_term_exp_ui(ulong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

    Set ``exp`` to the exponent vector of the term of index ``i``.

.. function:: void nmod_mpoly_set_term_exp_fmpz(nmod_mpoly_t A, slong i, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

.. function:: void nmod_mpoly_set_term_exp_ui(nmod_mpoly_t A, slong i, const ulong * exp, const nmod_mpoly_ctx_t ctx)

    Set the exponent of the term of index ``i`` to ``exp``.

.. function:: void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

.. function:: void nmod_mpoly_push_term_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)

    Append a term to ``A`` with coefficient ``c`` and exponent vector ``exp``.
    This function runs in constant average time.

.. function:: void nmod_mpoly_sort_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Sort the terms of ``A`` into the canonical ordering dictated by the ordering in ``ctx``.
    This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    This function runs in linear time in the bit size of ``A``.

.. function:: void nmod_mpoly_combine_like_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

    Combine adjacent like terms in ``A`` and delete terms with coefficient zero.
    If the terms of ``A`` were sorted to begin with, the result will be in canonical form.
    This function runs in linear time in the bit size of ``A``.

.. function:: void nmod_mpoly_reverse(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the reversal of ``B``.


Random generation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to  ``n_randint(state, exp_bound)``.

.. function:: void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bounds, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index ``i`` are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and exponents whose packed form does not exceed the given bit count.


Addition/Subtraction
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``c``.

.. function:: void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``c``.

.. function:: void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``C``.

.. function:: void nmod_mpoly_sub(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``C``.


Scalar operations
--------------------------------------------------------------------------------

.. function:: void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    
    Set ``A`` to `-```B``.

.. function:: void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``c``.

.. function:: void nmod_mpoly_make_monic(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` divided by the leading coefficient of ``B``.
    This throws if ``B`` is zero or the leading coefficient is not invertible.


Differentiation
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_derivative(nmod_mpoly_t A, const nmod_mpoly_t B, slong idx, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the derivative of ``B`` with respect to the variable of index ``idx``.


Evaluation
--------------------------------------------------------------------------------


.. function:: ulong nmod_mpoly_evaluate_all_ui(nmod_mpoly_t A, const ulong * vals, const nmod_mpoly_ctx_t ctx)

    Return the evaluation of ``A`` where the variables are replaced by the corresponding elements of the array ``vals``.

.. function:: void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong var, ulong val, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``var`` is replaced by ``val``.

.. function:: void nmod_mpoly_compose_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.

.. function:: void nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    Both ``A`` and the elements of ``C`` have context object ``ctxAC``, while ``B`` has context object ``ctxB``.
    Neither of ``A`` and ``B`` is allowed to alias any other polynomial.

    The compose functions try to guard against unreasonable arithmetic by throwing.


Multiplication
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_mul(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C``.

.. function:: void nmod_mpoly_mul_johnson(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C`` using Johnson's heap-based method.

.. function:: void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C`` using a heap and multiple threads.
    This function should only be called once ``global_thread_pool`` has been initialized.

.. function:: int nmod_mpoly_mul_array(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Try to set ``A`` to ``B`` times ``C`` using arrays.
    If the return is ``0``, the operation was unsuccessful. Otherwise, it was successful and the return is ``1``.

.. function:: int nmod_mpoly_mul_dense(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

    Try to set ``A`` to ``B`` times `C` using univariate arithmetic.
    If the return is ``0``, the operation was unsuccessful. Otherwise, it was successful and the return is ``1``.



Powering
--------------------------------------------------------------------------------


.. function:: void nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_t k, const nmod_mpoly_ctx_t ctx)

    Set `A` to `B` raised to the `k`-th power.
    This function throws if `k < 0` or if `k` does not fit an ``slong`` and `A` has more than one term.

.. function:: void nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)

    Set `A` to `B` raised to the `k`-th power.

.. function:: void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)

    Set `A` to `B` raised to the `k`-th power using repeated multiplications.


Division
--------------------------------------------------------------------------------

The division functions will generally throw if the leading coefficient of a divisor polynomial is not invertible.

.. function:: int nmod_mpoly_divides(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    If ``A`` is divisible by ``B``, set ``Q`` to the exact quotient and return ``1``. Otherwise, set ``Q`` to zero and return ``0``.
    Note that the function ``nmod_mpoly_div`` below may be faster if the quotient is known to be exact.

.. function:: void nmod_mpoly_div(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Set ``Q`` to the quotient of ``A`` by ``B``, discarding the remainder.

.. function:: void nmod_mpoly_divrem(nmod_mpoly_t Q, nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Set ``Q`` and ``R`` to the quotient and remainder of ``A`` divided by ``B``.

.. function:: void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R, const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len, const nmod_mpoly_ctx_t ctx)

    This function is as per :func:`fmpq_mpoly_divrem` except that it takes an array of divisor polynomials ``B`` and it returns an array of quotient polynomials ``Q``.
    The number of divisor (and hence quotient) polynomials, is given by ``len``.


.. function:: int nmod_mpoly_divides_dense(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Try to do the operation of ``nmod_mpoly_divides`` using univariate arithmetic.
    If the return is `-1`, the operation was unsuccessful. Otherwise, it was successful and the return is `0` or `1`.

.. function:: int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Do the operation of ``nmod_mpoly_divides`` using the algorithm of Michael Monagan and Roman Pearce.

.. function:: int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Do the operation of ``nmod_mpoly_divides`` using a heap and multiple threads.
    This function should only be called once ``global_thread_pool`` has been initialized.


Greatest Common Divisor
--------------------------------------------------------------------------------


.. function:: int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

    Try to set ``G`` to the monic GCD of ``A`` and ``B``. The GCD of zero and zero is defined to be zero.
    If the return is ``1`` the function was successful. Otherwise the return is  ``0`` and ``G`` is left untouched.
    If the modulus is not prime, this function will probably return ``0`` quickly.


Internal Functions
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


.. function:: void nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** q, nmod_mpoly_t r, const nmod_mpoly_t poly2, nmod_mpoly_struct * const * poly3, slong len, const nmod_mpoly_ctx_t ctx)

    This function is as per ``nmod_mpoly_divrem_monagan_pearce`` except
    that it takes an array of divisor polynomials ``poly3``, and it returns
    an array of quotient polynomials ``q``. The number of divisor (and hence
    quotient) polynomials, is given by ``len``. The function computes
    polynomials `q_i = q[i]` such that ``poly2`` is
    `r + \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where `b_i =` ``poly3[i]``.


.. function:: int nmod_mpoly_gcd_brown(nmod_mpoly_t poly1, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)

    If the return is nonzero, used Brown's dense modular algorithm to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` is monic.

.. function:: int nmod_mpoly_gcd_zippel(nmod_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

    If the return is nonzero, used a modular algorithm with Zippel's sparse
    interpolation to set
    ``poly1`` to the GCD of ``poly2`` and ``poly3``, where
    ``poly1`` is monic.

