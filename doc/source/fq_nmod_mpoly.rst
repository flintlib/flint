.. _fq_nmod-mpoly:

**fq_nmod_mpoly.h** -- multivariate polynomials over finite fields of word-sized characteristic
================================================================================================

    The exponents follow the ``mpoly`` interface.
    A coefficient may be referenced as a ``fq_nmod_struct *``.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_nmod_mpoly_ctx_struct

    Context structure for ``fq_nmod_mpoly``.

.. type:: fq_nmod_mpoly_ctx_t

    An array of length 1 of ``fq_nmod_mpoly_ctx_struct``.

.. type:: fq_nmod_mpoly_struct

    A structure holding a multivariate polynomial over a finite field of word-sized characteristic.

.. type:: fq_nmod_mpoly_t

    An array of length 1 of ``fq_nmod_mpoly_struct``.


Context object
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, const fq_nmod_ctx_t fqctx)

    Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    It will have coefficients in the finite field ``fqctx``.
    The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong fq_nmod_mpoly_ctx_nvars(fq_nmod_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: ordering_t fq_nmod_mpoly_ctx_ord(const fq_nmod_mpoly_ctx_t ctx)

    Return the ordering used to initialize the context.

.. function:: void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx)

    Release any space allocated by an ``fq_nmod_mpoly_ctx_t``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.

.. function:: void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc, const fq_nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

.. function:: void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and ``bits`` bits for the exponents.

.. function:: void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong len, const fq_nmod_mpoly_ctx_t ctx)

    Ensure that ``A`` has space for at least ``len`` terms.

.. function:: void fq_nmod_mpoly_fit_bits(fq_nmod_mpoly_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)

    Ensure that the exponent fields of ``A`` have at least ``bits`` bits.

.. function:: void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A, slong alloc, const fq_nmod_mpoly_ctx_t ctx)

    Reallocate ``A`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than ``alloc``.

.. function:: void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Release any space allocated for ``A``.


Input/Output
--------------------------------------------------------------------------------

    The variable strings in ``x`` start with the variable of most significance at index ``0``. If ``x`` is ``NULL``, the variables are named ``x1``, ``x2``, ect.

.. function:: char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing ``A``, given an array of variable strings ``x``.

.. function:: int fq_nmod_mpoly_fprint_pretty(FILE * file, const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``file``.

.. function:: int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``stdout``.

.. function:: int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the polynomial in the null-terminates string ``str`` given an array ``x`` of variable strings.
    If parsing ``str`` fails, ``A`` is set to zero, and ``-1`` is returned. Otherwise, ``0``  is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in ``x``. The character ``^`` must be immediately followed by the (integer) exponent.
    If any division is not exact, parsing fails.


Basic manipulation
--------------------------------------------------------------------------------

.. function:: void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the variable of index ``var``, where ``var = 0`` corresponds to the variable with the most significance with respect to the ordering. 

.. function:: int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

    If `var \ge 0`, return ``1`` if ``A`` is equal to the `var`-th generator, otherwise return ``0``.
    If `var < 0`, return ``1`` if the polynomial is equal to any generator, otherwise return ``0``.

.. function:: void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    
    Set ``A`` to ``B``.

.. function:: int fq_nmod_mpoly_equal(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to ``B``, else return ``0``.

.. function:: void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Efficiently swap ``A`` and ``B``.


Constants
--------------------------------------------------------------------------------


.. function:: int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is a constant, else return ``0``.

.. function:: void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Assuming that ``A`` is a constant, set ``c`` to this constant.
    This function throws if ``A`` is not a constant.

.. function:: void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``c``.

.. function:: void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant given by :func:`fq_nmod_gen`.

.. function:: void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``0``.

.. function:: void fq_nmod_mpoly_one(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``1``.

.. function:: int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``c``, else return ``0``.

.. function:: int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``0``, else return ``0``.

.. function:: int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``1``, else return ``0``.


Degrees
--------------------------------------------------------------------------------


.. function:: int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if the degrees of ``A`` with respect to each variable fit into an ``slong``, otherwise return ``0``.

.. function:: void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_degrees_si(slong * degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``A`` with respect to each variable.
    If ``A`` is zero, all degrees are set to ``-1``.

.. function:: void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
              slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Either return or set ``deg`` to the degree of ``A`` with respect to the variable of index ``var``.
    If ``A`` is zero, the degree is defined to be ``-1``.

.. function:: int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if the total degree of ``A`` fits into an ``slong``, otherwise return ``0``.

.. function:: void fq_nmod_mpoly_total_degree_fmpz(fmpz_t tdeg, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
              slong fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Either return or set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, the total degree is defined to be ``-1``.


Coefficients
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set ``c`` to the coefficient of the corresponding monomial in ``A``.
    This function thows if ``M`` is not a monomial.

.. function:: void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set the coefficient of the corresponding monomial in ``A`` to ``c``.
    This function thows if ``M`` is not a monomial.

.. function:: void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c, const fq_nmod_mpoly_t A, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c, const fq_nmod_mpoly_t A, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

    Set ``c`` to the coefficient of the monomial with exponent vector ``exp``.

.. function:: void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent ``exp`` to ``c``.

.. function:: void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C, const fq_nmod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fq_nmod_mpoly_ctx_t ctx)

    Set ``C`` to the coefficient of ``A`` with respect to the variables in ``vars`` with powers in the corresponding array ``exps``.
    Both ``vars`` and ``exps`` point to array of length ``length``. It is assumed that `0 < length \le nvars(A)` and that the variables in ``vars`` are distinct. 


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` (resp. ``-1``, or ``0``) if ``A`` is after (resp. before, same as) ``B`` in some arbitrary but fixed total ordering of the polynomials.
    This ordering agrees with the usual ordering of monomials when ``A`` and ``B`` are both monomials.


Container operations
--------------------------------------------------------------------------------

    These functions deal with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.

.. function:: fq_nmod_struct * fq_nmod_mpoly_term_coeff_ref(fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Return a reference to the coefficient of index `i` of ``A``.

.. function:: int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is in canonical form. Otherwise, return ``0``.
    To be in canonical form, all of the terms must have nonzero coefficients, and the terms must be sorted from greatest to least.

.. function:: slong fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return the number of terms in ``A``.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, slong new_length, const fq_nmod_mpoly_ctx_t ctx)

    Set the length of ``A`` to ``new_length``.
    Terms are either deleted from the end, or new zero terms are appended.

.. function:: void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set ``c`` to the coefficient of the term of index ``i``.

.. function:: void fq_nmod_mpoly_set_term_coeff_ui(fq_nmod_mpoly_t A, slong i, ulong c, const fq_nmod_mpoly_ctx_t ctx)

    Set the coefficient of the term of index ``i`` to ``c``.

.. function:: int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
              int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Return ``1`` if all entries of the exponent vector of the term of index `i` fit into an ``slong`` (resp. a ``ulong). Otherwise, return ``0``.

.. function:: void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_get_term_exp_ui(ulong * exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_get_term_exp_si(slong * exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set ``exp`` to the exponent vector of the term of index ``i``.
    The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

.. function:: ulong fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A, slong i, slong var, const fq_nmod_mpoly_ctx_t ctx)
              slong fq_nmod_mpoly_get_term_var_exp_si(const fq_nmod_mpoly_t A, slong i, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Return the exponent of the variable ``var`` of the term of index ``i``.
    This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

.. function:: void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A, slong i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A, slong i, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

    Set the exponent of the term of index ``i`` to ``exp``.

.. function:: void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set ``M`` to the term of index ``i`` in ``A``.

.. function:: void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set ``M`` to the monomial of the term of index ``i`` in ``A``. The coefficient of ``M`` will be one.

.. function:: void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

    Append a term to ``A`` with coefficient ``c`` and exponent vector ``exp``.
    This function runs in constant average time.

.. function:: void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Sort the terms of ``A`` into the canonical ordering dictated by the ordering in ``ctx``.
    This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    This function runs in linear time in the bit size of ``A``.

.. function:: void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

    Combine adjacent like terms in ``A`` and delete terms with coefficient zero.
    If the terms of ``A`` were sorted to begin with, the result will be in canonical form.
    This function runs in linear time in the bit size of ``A``.

.. function:: void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the reversal of ``B``.


Random generation
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const fq_nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to  ``n_randint(state, exp_bound)``.

.. function:: void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bounds, const fq_nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index ``i`` are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const fq_nmod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and exponents whose packed form does not exceed the given bit count.


Addition/Subtraction
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``c``.

.. function:: void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``c``.

.. function:: void fq_nmod_mpoly_add(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``C``.

.. function:: void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``C``.


Scalar operations
--------------------------------------------------------------------------------

.. function:: void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    
    Set ``A`` to `-```B``.

.. function:: void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``c``.

.. function:: void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` divided by the leading coefficient of ``B``.
    This throws if ``B`` is zero.


Differentiation
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the derivative of ``B`` with respect to the variable of index ``idx``.


Evaluation
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, fq_nmod_mpoly_t A, fq_nmod_struct * const *  vals, const fq_nmod_mpoly_ctx_t ctx)

    Set ``ev`` the evaluation of ``A`` where the variables are replaced by the corresponding elements of the array ``vals``.

.. function:: void fq_nmod_mpoly_evaluate_one_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, slong var, fq_nmod_t val, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``var`` is replaced by ``val``.

.. function:: int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A, const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.
    Return `1` for success and `0` for failure.

.. function:: int fq_nmod_mpoly_compose_fq_nmod_mpoly(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    Both ``A`` and the elements of ``C`` have context object ``ctxAC``, while ``B`` has context object ``ctxB``.
    Neither ``A`` nor ``B`` is allowed to alias any other polynomial.
    Return `1` for success and `0` for failure.

.. function:: void fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const slong * c, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``i`` in ``ctxB`` is replaced by the variable of index ``c[i]`` in ``ctxAC``.
    The length of the array ``C`` is the number of variables in ``ctxB``.
    If any ``c[i]`` is negative, the corresponding variable of ``B`` is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in ``ctxAC``.


Multiplication
--------------------------------------------------------------------------------


.. function:: void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C``.


Powering
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx)

    Set `A` to `B` raised to the `k`-th power.
    Return `1` for success and `0` for failure.

.. function:: int fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, ulong k, const fq_nmod_mpoly_ctx_t ctx)

    Set `A` to `B` raised to the `k`-th power.
    Return `1` for success and `0` for failure.


Division
--------------------------------------------------------------------------------


.. function:: int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    If ``A`` is divisible by ``B``, set ``Q`` to the exact quotient and return ``1``. Otherwise, set ``Q`` to zero and return ``0``.

.. function:: void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Set ``Q`` to the quotient of ``A`` by ``B``, discarding the remainder.

.. function:: void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Set ``Q`` and ``R`` to the quotient and remainder of ``A`` divided by ``B``.

.. function:: void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, fq_nmod_mpoly_struct * const * B, slong len, const fq_nmod_mpoly_ctx_t ctx)

    This function is as per :func:`fq_nmod_mpoly_divrem` except that it takes an array of divisor polynomials ``B`` and it returns an array of quotient polynomials ``Q``.
    The number of divisor (and hence quotient) polynomials, is given by ``len``.


Greatest Common Divisor
--------------------------------------------------------------------------------


.. function:: int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

    Try to set ``G`` to the monic GCD of ``A`` and ``B``. The GCD of zero and zero is defined to be zero.
    If the return is ``1`` the function was successful. Otherwise the return is  ``0`` and ``G`` is left untouched.


Univariate Functions
--------------------------------------------------------------------------------

    An ``fq_nmod_mpoly_univar_t`` holds a univariate polynomial in some main variable
    with ``fq_nmod_mpoly_t`` coefficients in the remaining variables. These functions
    are useful when one wants to rewrite an element of `\mathbb{F}_q[x_1, \dots, x_m]`
    as an element of `(\mathbb{F}_q[x_1, \dots, x_{v-1}, x_{v+1}, \dots, x_m])[x_v]`
    and vise versa.

.. function:: void fq_nmod_mpoly_univar_init(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)

    Initialize `A`.

.. function:: void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)

    Clear `A`.

.. function:: void fq_nmod_mpoly_univar_swap(fq_nmod_mpoly_univar_t A, fq_nmod_mpoly_univar_t B, const fq_nmod_mpoly_ctx_t ctx)

    Swap `A` and `B`.

.. function:: void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to a univariate form of ``B`` by pulling out the variable of index ``var``.
    The coefficients of ``A`` will still belong to the content ``ctx`` but will not depend on the variable of index ``var``.

.. function:: void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

    Set ``A`` to the normal form of ``B`` by putting in the variable of index ``var``.
    This function is undefined if the coefficients of ``B`` depend on the variable of index ``var``.

.. function:: int fq_nmod_mpoly_univar_degree_fits_si(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return `1` if the degree of ``A`` with respect to the main variable fits an ``slong``. Otherwise, return `0`.

.. function:: slong fq_nmod_mpoly_univar_length(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)

    Return the number of terms in ``A`` with respect to the main variable.

.. function:: slong fq_nmod_mpoly_univar_get_term_exp_si(fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Return the exponent of the term of index ``i`` of ``A``.

.. function:: void fq_nmod_mpoly_univar_get_term_coeff(fq_nmod_mpoly_t c, const fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
              void fq_nmod_mpoly_univar_swap_term_coeff(fq_nmod_mpoly_t c, fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

    Set (resp. swap) ``c`` to (resp. with) the coefficient of the term of index ``i`` of ``A``.


