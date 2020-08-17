.. _fmpq-mpoly:

**fmpq_mpoly.h** -- multivariate polynomials over the rational numbers
===============================================================================

    The exponents follow the ``mpoly`` interface.
    No references to the coefficients are available.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpq_mpoly_ctx_struct

    Context structure for ``fmpq_mpoly``.

.. type:: fmpq_mpoly_ctx_t

    An array of length 1 of ``fmpq_mpoly_ctx_struct``.

.. type:: fmpq_mpoly_struct

    A structure holding a multivariate rational polynomial.  It is implemented as a
    ``fmpq_t`` holding the content of the polynomial and a primitive integer
    polynomial.

.. type:: fmpq_mpoly_t

    An array of length 1 of ``fmpq_mpoly_struct``.


Context object
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

    Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong fmpq_mpoly_ctx_nvars(const fmpq_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: ordering_t fmpq_mpoly_ctx_ord(const fmpq_mpoly_ctx_t ctx)

    Return the ordering used to initialize the context.

.. function:: void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)

    Release up any space allocated by an ``fmpq_mpoly_ctx_t``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_init(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.

.. function:: void fmpq_mpoly_init2(fmpq_mpoly_t A, slong alloc, const fmpq_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

.. function:: void fmpq_mpoly_init3(fmpq_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and ``bits`` bits for the exponents.

.. function:: void fmpq_mpoly_fit_length(fmpq_mpoly_t A, slong len, const fmpq_mpoly_ctx_t ctx)

    Ensure that ``A`` has space for at least ``len`` terms.

.. function:: void fmpq_mpoly_fit_bits(fmpq_mpoly_t A, flint_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)

    Ensure that the exponent fields of ``A`` have at least ``bits`` bits.

.. function:: void fmpq_mpoly_realloc(fmpq_mpoly_t A, slong alloc, const fmpq_mpoly_ctx_t ctx)

    Reallocate ``A`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than ``alloc``.

.. function:: void fmpq_mpoly_clear(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Release any space allocated for ``A``.


Input/Output
--------------------------------------------------------------------------------

    The variable strings in ``x`` start with the variable of most significance at index ``0``. If ``x`` is ``NULL``, the variables are named ``x1``, ``x2``, ect.

.. function:: char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing ``A``, given an array of variable strings ``x``.

.. function:: int fmpq_mpoly_fprint_pretty(FILE * file, const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``file``.

.. function:: int fmpq_mpoly_print_pretty(const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``stdout``.

.. function:: int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char * str, const char ** x, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the polynomial in the null-terminates string ``str`` given an array ``x`` of variable strings.
    If parsing ``str`` fails, ``A`` is set to zero, and ``-1`` is returned. Otherwise, ``0``  is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in ``x``. The character ``^`` must be immediately followed by the (integer) exponent.
    If any division is not exact, parsing fails.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_gen(fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the variable of index ``var``, where ``var = 0`` corresponds to the variable with the most significance with respect to the ordering. 

.. function:: int fmpq_mpoly_is_gen(const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

    If `var \ge 0`, return ``1`` if ``A`` is equal to the `var`-th generator, otherwise return ``0``.
    If `var < 0`, return ``1`` if the polynomial is equal to any generator, otherwise return ``0``.

.. function:: void fmpq_mpoly_set(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
    
    Set ``A`` to ``B``.

.. function:: int fmpq_mpoly_equal(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to ``B``, else return ``0``.

.. function:: void fmpq_mpoly_swap(fmpq_mpoly_t A, fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Efficiently swap ``A`` and ``B``.


Constants
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is a constant, else return ``0``.

.. function:: void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``A`` is a constant, set ``c`` to this constant.
    This function throws if ``A`` is not a constant.

.. function:: void fmpq_mpoly_set_fmpq(fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_set_fmpz(fmpq_mpoly_t A, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_set_ui(fmpq_mpoly_t A, ulong c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_set_si(fmpq_mpoly_t A, slong c, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``c``.

.. function:: void fmpq_mpoly_zero(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``0``.

.. function:: void fmpq_mpoly_one(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``1``.

.. function:: int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t A, fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t A, fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_equal_ui(const fmpq_mpoly_t A, ulong c, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_equal_si(const fmpq_mpoly_t A, slong c, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``c``, else return ``0``.

.. function:: int fmpq_mpoly_is_zero(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``0``, else return ``0``.

.. function:: int fmpq_mpoly_is_one(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``1``, else return ``0``.


Degrees
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if the degrees of ``A`` with respect to each variable fit into an ``slong``, otherwise return ``0``.

.. function:: void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``degs`` to the degrees of ``A`` with respect to each variable.
    If ``A`` is zero, all degrees are set to ``-1``.

.. function:: void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)
              slong fmpq_mpoly_degree_si(const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

    Either return or set ``deg`` to the degree of ``A`` with respect to the variable of index ``var``.
    If ``A`` is zero, the degree is defined to be ``-1``.

.. function:: int fmpq_mpoly_total_degree_fits_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if the total degree of ``A`` fits into an ``slong``, otherwise return ``0``.

.. function:: void fmpq_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
              slong fmpq_mpoly_total_degree_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Either return or set ``tdeg`` to the total degree of ``A``.
    If ``A`` is zero, the total degree is defined to be ``-1``.


Coefficients
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_get_denominator(fmpz_t d, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``d`` to the denominator of ``A``, the smallest positive integer `d` such that `d*A` has integer coefficients.

.. function:: void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c, const fmpq_mpoly_t A, const fmpq_mpoly_t M, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set ``c`` to the coefficient of the corresponding monomial in ``A``.
    This function thows if ``M`` is not a monomial.

.. function:: void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_t M, const fmpq_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set the coefficient of the corresponding monomial in ``A`` to ``c``.
    This function thows if ``M`` is not a monomial.

.. function:: void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t A, ulong const * exp, const fmpq_mpoly_ctx_t ctx)

    Set ``c`` to the coefficient of the monomial with exponent ``exp``.

.. function:: void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A, const fmpq_t c, fmpz * const * exp, fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t A, const fmpq_t c, ulong const * exp, fmpq_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent ``exp`` to ``c``.

.. function:: void fmpq_mpoly_get_coeff_vars_ui(fmpq_mpoly_t C, const fmpq_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpq_mpoly_ctx_t ctx)

    Set ``C`` to the coefficient of ``A`` with respect to the variables in ``vars`` with powers in the corresponding array ``exps``.
    Both ``vars`` and ``exps`` point to array of length ``length``. It is assumed that `0 < length \le nvars(A)` and that the variables in ``vars`` are distinct. 


Comparison
--------------------------------------------------------------------------------

.. function:: int fmpq_mpoly_cmp(const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` (resp. ``-1``, or ``0``) if ``A`` is after (resp. before, same as) ``B`` in some arbitrary but fixed total ordering of the polynomials.
    This ordering agrees with the usual ordering of monomials when ``A`` and ``B`` are both monomials.


Container operations
--------------------------------------------------------------------------------

    These function try to deal efficiently with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.
    The mutating functions here are not guaranteed to leave the polynomial in reduced form (see :func:`fmpq_mpoly_is_canonical` for a definition of reduced).
    This means that even if nonzero terms with distinct exponents have been constructed in the correct order, a call to :func:`fmpq_mpoly_reduce` is necessary to ensure that the polynomial is in canonical form.
    As with the ``fmpz_mpoly`` module, a call to :func:`fmpq_mpoly_sort_terms` followed by a call to :func:`fmpq_mpoly_combine_like_terms` should leave the polynomial in canonical form.

.. function:: fmpq * fmpq_mpoly_content_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return a reference to the content of ``A``.

.. function:: fmpz_mpoly_struct * fmpq_mpoly_zpoly_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return a reference to the integer polynomial of ``A``.

.. function:: fmpz * fmpq_mpoly_zpoly_term_coeff_ref(fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Return a reference to the coefficient of index `i` of the integer polynomial of ``A``.

.. function:: int fmpq_mpoly_is_canonical(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is in canonical form. Otherwise, return ``0``.
    An ``fmpq_mpoly_t`` is represented as the product of an ``fmpq_t content`` and an ``fmpz_mpoly_t zpoly``.
    The representation is considered canonical when either
    (1) both ``content`` and ``zpoly`` are zero, or
    (2) both ``content`` and ``zpoly`` are nonzero and canonical and ``zpoly`` is reduced.
    A nonzero ``zpoly`` is considered reduced when the coefficients have GCD one and the leading coefficient is positive.

.. function:: slong fmpq_mpoly_length(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Return the number of terms stored in ``A``.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void fmpq_mpoly_resize(fmpq_mpoly_t A, slong new_length, const fmpq_mpoly_ctx_t ctx)

    Set the length of ``A`` to ``new_length``.
    Terms are either deleted from the end, or new zero terms are appended.

.. function:: void fmpq_mpoly_get_term_coeff_fmpq(fmpq_t c, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Set `c` to coefficient of index `i`

.. function:: void fmpq_mpoly_set_term_coeff_fmpq(fmpq_mpoly_t A, slong i, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)

    Set the coefficient of index `i` to `c`.

.. function:: int fmpq_mpoly_term_exp_fits_si(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_term_exp_fits_ui(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Return ``1`` if all entries of the exponent vector of the term of index `i`  fit into an ``slong`` (resp. a ``ulong``). Otherwise, return ``0``.

.. function:: void fmpq_mpoly_get_term_exp_fmpz(fmpz ** exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_get_term_exp_ui(ulong * exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_get_term_exp_si(slong * exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Set ``exp`` to the exponent vector of the term of index ``i``.
    The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

.. function:: ulong fmpq_mpoly_get_term_var_exp_ui(const fmpq_mpoly_t A, slong i, slong var, const fmpq_mpoly_ctx_t ctx)
              slong fmpq_mpoly_get_term_var_exp_si(const fmpq_mpoly_t A, slong i, slong var, const fmpq_mpoly_ctx_t ctx)

    Return the exponent of the variable ``var`` of the term of index ``i``.
    This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

.. function:: void fmpq_mpoly_set_term_exp_fmpz(fmpq_mpoly_t A, slong i, fmpz * const * exps, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_set_term_exp_ui(fmpq_mpoly_t A, slong i, const ulong * exps, const fmpq_mpoly_ctx_t ctx)

    Set the exponent vector of the term of index ``i`` to ``exp``.

.. function:: void fmpq_mpoly_get_term(fmpq_mpoly_t M, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Set ``M`` to the term of index ``i`` in ``A``.

.. function:: void fmpq_mpoly_get_term_monomial(fmpq_mpoly_t M, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Set ``M`` to the monomial of the term of index ``i`` in ``A``. The coefficient of ``M`` will be one.

.. function:: void fmpq_mpoly_push_term_fmpq_fmpz(fmpq_mpoly_t A, const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_fmpz_fmpz(fmpq_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_ui_fmpz(fmpq_mpoly_t A, ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_si_fmpz(fmpq_mpoly_t A, slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_t A, const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_fmpz_ui(fmpq_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_ui_ui(fmpq_mpoly_t A, ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_push_term_si_ui(fmpq_mpoly_t A, slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

    Append a term to ``A`` with coefficient ``c`` and exponent vector ``exp``.
    This function should run in constant average time if the terms pushed have bounded denominator.

.. function:: void fmpq_mpoly_reduce(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Factor out necessary content from ``A->zpoly`` so that it is reduced.
    If the terms of ``A`` were nonzero and sorted with distinct exponents to begin with, the result will be in canonical form.

.. function:: void fmpq_mpoly_sort_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Sort the internal ``A->zpoly`` into the canonical ordering dictated by the ordering in ``ctx``.
    This function does not combine like terms, nor does it delete terms with coefficient zero, nor does it reduce.

.. function:: void fmpq_mpoly_combine_like_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Combine adjacent like terms in the internal ``A->zpoly`` and then factor out content via a call to :func:`fmpq_mpoly_reduce`.
    If the terms of ``A`` were sorted to begin with, the result will be in canonical form.

.. function:: void fmpq_mpoly_reverse(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the reversal of ``B``.


Random generation
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_randtest_bound(fmpq_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong exp_bound, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to ``n_randint(state, exp_bound)``.

.. function:: void fmpq_mpoly_randtest_bounds(fmpq_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong * exp_bounds, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index ``i`` are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void fmpq_mpoly_randtest_bits(fmpq_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, mp_limb_t exp_bits, const fmpq_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and exponents whose packed form does not exceed the given bit count.

    The parameter ``coeff_bits`` to the three functions ``fmpq_mpoly_randtest_{bound|bounds|bits}`` is merely a suggestion for the approximate bit count of the resulting coefficients.


Addition/Subtraction
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_add_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_add_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_add_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_add_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``c``.

.. function:: void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_sub_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_sub_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``c``.

.. function:: void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``C``.

.. function:: void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``C``.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_neg(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
    
    Set ``A`` to `-```B``.

.. function:: void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times `c`.

.. function:: void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` divided by `c`.

.. function:: void fmpq_mpoly_make_monic(fmpq_mpoly_t A, fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` divided by the leading coefficient of ``B``.
    This throws if ``B`` is zero.

    All of these functions run quickly if ``A`` and ``B`` are aliased.


Differentiation/Integration
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_derivative(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the derivative of ``B`` with respect to the  variable of index ``var``.

.. function:: void fmpq_mpoly_integral(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the integral with the fewest number of terms of ``B`` with respect to the variable of index ``var``.


Evaluation
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A, fmpq * const * vals, const fmpq_mpoly_ctx_t ctx)

    Set ``ev`` the evaluation of ``A`` where the variables are replaced by the corresponding elements of the array ``vals``.
    Return `1` for success and `0` for failure.

.. function:: int fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_t val, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``var`` is replaced by ``val``.
    Return `1` for success and `0` for failure.

.. function:: int fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A, const fmpq_mpoly_t B, fmpq_poly_struct * const * C, const fmpq_mpoly_ctx_t ctxB)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.
    Return `1` for success and `0` for failure.

.. function:: int fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A, const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C, const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    Both ``A`` and the elements of ``C`` have context object ``ctxAC``, while ``B`` has context object ``ctxB``.
    Neither ``A`` nor ``B`` is allowed to alias any other polynomial.
    Return `1` for success and `0` for failure.

.. function:: void fmpq_mpoly_compose_fmpq_mpoly_gen(fmpq_mpoly_t A, const fmpq_mpoly_t B, const slong * c, const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``i`` in ``ctxB`` is replaced by the variable of index ``c[i]`` in ``ctxAC``.
    The length of the array ``C`` is the number of variables in ``ctxB``.
    If any ``c[i]`` is negative, the corresponding variable of ``B`` is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in ``ctxAC``.


Multiplication
--------------------------------------------------------------------------------


.. function:: void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C``.


Powering
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t k, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` raised to the `k`-th power.
    Return `1` for success and `0` for failure.

.. function:: int fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong k, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` raised to the `k`-th power.
    Return `1` for success and `0` for failure.


Division
--------------------------------------------------------------------------------


.. function:: int fmpq_mpoly_divides(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    If ``A`` is divisible by ``B``, set ``Q`` to the exact quotient and return ``1``. Otherwise, set ``Q`` to zero and return ``0``.
    Note that the function :func:`fmpq_mpoly_div` may be faster if the quotient is known to be exact.

.. function:: void fmpq_mpoly_div(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Set ``Q`` to the quotient of ``A`` by ``B``, discarding the remainder.

.. function:: void fmpq_mpoly_divrem(fmpq_mpoly_t Q, fmpq_mpoly_t R, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Set ``Q`` and ``R`` to the quotient and remainder of ``A`` divided by ``B``.

.. function:: void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** Q, fmpq_mpoly_t R, const fmpq_mpoly_t A, fmpq_mpoly_struct * const * B, slong len, const fmpq_mpoly_ctx_t ctx)

    This function is as per :func:`fmpq_mpoly_divrem` except that it takes an array of divisor polynomials ``B`` and it returns an array of quotient polynomials ``Q``.
    The number of divisor (and hence quotient) polynomials, is given by ``len``.


Greatest Common Divisor
--------------------------------------------------------------------------------

.. function:: void fmpq_mpoly_content(fmpq_t g, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``g`` to the (nonnegative) gcd of the coefficients of ``A``.

.. function:: void fmpq_mpoly_term_content(fmpq_mpoly_t M, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

    Set ``M`` to the GCD of the terms of ``A``.
    If ``A`` is zero, ``M`` will be zero. Otherwise, ``M`` will be a monomial with coefficient one.

.. function:: int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
              int fmpq_mpoly_gcd_threaded(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx, slong thread_limit)

    Try to set ``G`` to the monic GCD of ``A`` and ``B``. The GCD of zero and zero is defined to be zero.
    If the return is ``1`` the function was successful. Otherwise the return is  ``0`` and ``G`` is left untouched.

.. function:: int fmpq_mpoly_gcd_cofactors(fmpq_mpoly_t G, fmpq_mpoly_t Abar, fmpq_mpoly_t Bbar, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

    Do the operation of :func:`fmpq_mpoly_gcd` and also compute ``Abar = A/G`` and ``Bbar = B/G`` if successful.


Univariate Functions
--------------------------------------------------------------------------------

    An ``fmpq_mpoly_univar_t`` holds a univariate polynomial in some main variable
    with ``fmpq_mpoly_t`` coefficients in the remaining variables. These functions
    are useful when one wants to rewrite an element of `\mathbb{Q}[x_1, \dots, x_m]`
    as an element of `(\mathbb{Q}[x_1, \dots, x_{v-1}, x_{v+1}, \dots, x_m])[x_v]`
    and vise versa.

.. function:: void fmpq_mpoly_univar_init(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)

    Initialize `A`.

.. function:: void fmpq_mpoly_univar_clear(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)

    Clear `A`.

.. function:: void fmpq_mpoly_univar_swap(fmpq_mpoly_univar_t A, fmpq_mpoly_univar_t B, const fmpq_mpoly_ctx_t ctx)

    Swap `A` and `B`.

.. function:: void fmpq_mpoly_to_univar(fmpq_mpoly_univar_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to a univariate form of ``B`` by pulling out the variable of index ``var``.
    The coefficients of ``A`` will still belong to the content ``ctx`` but will not depend on the variable of index ``var``.

.. function:: void fmpq_mpoly_from_univar(fmpq_mpoly_t A, const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx)

    Set ``A`` to the normal form of ``B`` by putting in the variable of index ``var``.
    This function is undefined if the coefficients of ``B`` depend on the variable of index ``var``.

.. function:: int fmpq_mpoly_univar_degree_fits_si(const fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)

    Return `1` if the degree of ``A`` with respect to the main variable fits an ``slong``. Otherwise, return `0`.

.. function:: slong fmpq_mpoly_univar_length(const fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)

    Return the number of terms in ``A`` with respect to the main variable.

.. function:: slong fmpq_mpoly_univar_get_term_exp_si(fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Return the exponent of the term of index ``i`` of ``A``.

.. function:: void fmpq_mpoly_univar_get_term_coeff(fmpq_mpoly_t c, const fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)
              void fmpq_mpoly_univar_swap_term_coeff(fmpq_mpoly_t c, fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)

    Set (resp. swap) ``c`` to (resp. with) the coefficient of the term of index ``i`` of ``A``.


