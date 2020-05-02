.. _fmpz-mpoly:

**fmpz_mpoly.h** -- multivariate polynomials over the integers
===============================================================================

    The exponents follow the ``mpoly`` interface.
    A coefficient may be referenced as a ``fmpz *``.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mpoly_struct

    Context structure for ``fmpz_mpoly``.

.. type:: fmpz_mpoly_t

    An array of length 1 of ``fmpz_mpoly_ctx_struct``.

.. type:: fmpz_mpoly_ctx_struct

    A structure holding a multivariate integer polynomial.

.. type:: fmpz_mpoly_ctx_t

    An array of length 1 of ``fmpz_mpoly_struct``.


Context object
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

    Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong fmpz_mpoly_ctx_nvars(fmpz_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: ordering_t fmpz_mpoly_ctx_ord(const fmpz_mpoly_ctx_t ctx)

    Return the ordering used to initialize the context.

.. function:: void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)

    Release up any space allocated by an ``fmpz_mpoly_ctx_t``.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.

.. function:: void fmpz_mpoly_init2(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

.. function:: void fmpz_mpoly_init3(fmpz_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)

    Initialise ``A`` for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for ``alloc`` terms and ``bits`` bits for the exponents.

.. function:: void fmpz_mpoly_fit_length(fmpz_mpoly_t A, slong len, const fmpz_mpoly_ctx_t ctx)

    Ensure that ``A`` has space for at least ``len`` terms.

.. function:: void fmpz_mpoly_fit_bits(fmpz_mpoly_t A, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)

    Ensure that the exponent fields of ``A`` have at least ``bits`` bits.

.. function:: void fmpz_mpoly_realloc(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx)

    Reallocate ``A`` to have space for ``alloc`` terms. 
    Assumes the current length of the polynomial is not greater than ``alloc``.

.. function:: void fmpz_mpoly_clear(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Release any space allocated for ``A``.


Input/Output
--------------------------------------------------------------------------------

    The variable strings in ``x`` start with the variable of most significance at index ``0``. If ``x`` is ``NULL``, the variables are named ``x1``, ``x2``, ect.

.. function:: char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing ``A``, given an array of variable strings ``x``.

.. function:: int fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``file``.

.. function:: int fmpz_mpoly_print_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Print a string representing ``A`` to ``stdout``.

.. function:: int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t A, const char * str, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the polynomial in the null-terminates string ``str`` given an array ``x`` of variable strings.
    If parsing ``str`` fails, ``A`` is set to zero, and ``-1`` is returned. Otherwise, ``0``  is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in ``x``. The character ``^`` must be immediately followed by the (integer) exponent.
    If any division is not exact, parsing fails.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_gen(fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the variable of index ``var``, where ``var = 0`` corresponds to the variable with the most significance with respect to the ordering. 

.. function:: int fmpz_mpoly_is_gen(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

    If `var \ge 0`, return ``1`` if ``A`` is equal to the `var`-th generator, otherwise return ``0``.
    If `var < 0`, return ``1`` if the polynomial is equal to any generator, otherwise return ``0``.

.. function:: void fmpz_mpoly_set(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    
    Set ``A`` to ``B``.

.. function:: int fmpz_mpoly_equal(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to ``B``, else return ``0``.

.. function:: void fmpz_mpoly_swap(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Efficiently swap ``A`` and ``B``.

.. function:: int _fmpz_mpoly_fits_small(const fmpz * poly, slong len)

    Return 1 if the array of coefficients of length ``len`` consists
    entirely of values that are small ``fmpz`` values, i.e. of at most
    ``FLINT_BITS - 2`` bits plus a sign bit.

.. function:: slong fmpz_mpoly_max_bits(const fmpz_mpoly_t A)

    Computes the maximum number of bits `b` required to represent the absolute
    values of the coefficients of ``A``. If all of the coefficients are
    positive, `b` is returned, otherwise `-b` is returned.


Constants
--------------------------------------------------------------------------------


.. function:: int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is a constant, else return ``0``.

.. function:: void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``A`` is a constant, set ``c`` to this constant.
    This function throws if ``A`` is not a constant.

.. function:: void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_ui(fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``c``.

.. function:: void fmpz_mpoly_zero(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``0``.

.. function:: void fmpz_mpoly_one(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the constant ``1``.

.. function:: int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t A, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: int fmpz_mpoly_equal_si(const fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is equal to the constant ``c``, else return ``0``.

.. function:: int fmpz_mpoly_is_zero(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``0``, else return ``0``.

.. function:: int fmpz_mpoly_is_one(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is the constant ``1``, else return ``0``.


Degrees
--------------------------------------------------------------------------------


.. function:: int fmpz_mpoly_degrees_fit_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

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
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set ``c`` to the coefficient of the corresponding monomial in ``A``.
    This function thows if ``M`` is not a monomial.

.. function:: void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t poly, const fmpz_t c, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

    Assuming that ``M`` is a monomial, set the coefficient of the corresponding monomial in ``A`` to ``c``.
    This function thows if ``M`` is not a monomial.

.. function:: void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: ulong fmpz_mpoly_get_coeff_ui_fmpz(const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_get_coeff_si_fmpz(const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: ulong fmpz_mpoly_get_coeff_ui_ui(const fmpz_mpoly_t A, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_get_coeff_si_ui(const fmpz_mpoly_t A, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

    Either return or set ``c`` to the coefficient of the monomial with exponent vector ``exp``.

.. function:: void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, ulong const * exp, fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t A, ulong c, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t A, slong c, ulong const * exp, const fmpz_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector ``exp`` to ``c``.

.. function:: void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C, const fmpz_mpoly_t A, slong * vars, ulong * exps, slong length, const fmpz_mpoly_ctx_t ctx)

    Set ``C`` to the coefficient of ``A`` with respect to the variables in ``vars`` with powers in the corresponding array ``exps``.
    Both ``vars`` and ``exps`` point to array of length ``length``. It is assumed that `0 < length \le nvars(A)` and that the variables in ``vars`` are distinct.


Comparison
--------------------------------------------------------------------------------


.. function:: int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` (resp. ``-1``, or ``0``) if the monomial of ``A`` is greater than (resp. less than, same as) the monomial of ``B``.
    ``A`` and ``B`` should both have length one with coefficient one. This function will throw otherwise.


Container operations
--------------------------------------------------------------------------------

    These functions deal with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.

.. function:: fmpz * fmpz_mpoly_term_coeff_ref(fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Return a reference to the coefficient of index `i` of ``A``.

.. function:: int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return ``1`` if ``A`` is in canonical form. Otherwise, return ``0``.
    To be in canonical form, all of the terms must have nonzero coefficient, and the terms must be sorted from greatest to least.

.. function:: slong fmpz_mpoly_length(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Return the number of terms in ``A``.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void fmpz_mpoly_resize(fmpz_mpoly_t A, slong new_length, const fmpz_mpoly_ctx_t ctx)

    Set the length of ``A`` to ``new_length``.
    Terms are either deleted from the end, or new zero terms are appended.

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

.. function:: void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Set ``exp`` to the exponent vector of the term of index ``i``.
    The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

.. function:: ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx)

.. function:: slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx)

    Return the exponent of the variable ``var`` of the term of index ``i``.
    This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

.. function:: void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A, slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    Set the exponent vector of the term of index ``i`` to ``exp``.

.. function:: void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Set ``M`` to the term of index ``i`` in ``A``.

.. function:: void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Set ``M`` to the monomial of the term of index ``i`` in ``A``. The coefficient of ``M`` will be one.

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
    This function runs in linear time in the size of ``A``.

.. function:: void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Combine adjacent like terms in ``A`` and delete terms with coefficient zero.
    If the terms of ``A`` were sorted to begin with, the result will be in canonical form.
    This function runs in linear time in the size of ``A``.

.. function:: void fmpz_mpoly_reverse(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the reversal of ``B``.


Random generation
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_randtest_bound(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong exp_bound, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to ``n_randint(state, exp_bound)``.

.. function:: void fmpz_mpoly_randtest_bounds(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, ulong * exp_bounds, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to ``length`` and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index ``i`` are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state, slong length, mp_limb_t coeff_bits, mp_limb_t exp_bits, const fmpz_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to the given length and exponents whose packed form does not exceed the given bit count.

    The parameter ``coeff_bits`` to the three functions ``fmpz_mpoly_randtest_{bound|bounds|bits}`` is merely a suggestion for the approximate bit count of the resulting signed coefficients.
    The function :func:`fmpz_mpoly_max_bits` will give the exact bit count of the result.


Addition/Subtraction
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_add_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus `c`.
    If ``A`` and ``B`` are aliased, this function will probably run quickly.

.. function:: void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_sub_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_sub_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus `c`.
    If ``A`` and ``B`` are aliased, this function will probably run quickly.

.. function:: void fmpz_mpoly_add(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` plus ``C``.
    If ``A`` and ``B`` are aliased, this function might run in time proportional to the size of ``C``.
    
.. function:: void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` minus ``C``.
    If ``A`` and ``B`` are aliased, this function might run in time proportional to the size of ``C``.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B const fmpz_mpoly_ctx_t ctx)
    
    Set ``A`` to ``-B``.

.. function:: void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``c``.

.. function:: void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` divided by ``c``. The division is assumed to be exact.

.. function:: int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)

.. function:: int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

.. function:: int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

    If ``B`` is divisible by ``c``, set ``A`` to the exact quotient and return ``1``, otherwise set ``A`` to zero and return ``0``.


Differentiation/Integration
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_derivative(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the derivative of ``B`` with respect to the variable of index ``var``.

.. function:: void fmpz_mpoly_integral(fmpz_mpoly_t A, fmpz_t scale, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` and ``scale`` so that ``A`` is an integral of ``scale*B`` with respect to the variable of index ``idx``, where ``scale`` is positive and as small as possible.


Evaluation
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A, fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)

    Set ``ev`` to the evaluation of ``A`` where the variables are replaced by the corresponding elements of the array ``vals``.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``var`` is replaced by ``val``.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctxB)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    The context object of ``B`` is ``ctxB``.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mpoly_compose_fmpz_mpoly_geobucket(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

.. function:: int fmpz_mpoly_compose_fmpz_mpoly_horner(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

.. function:: int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variables are replaced by the corresponding elements of the array ``C``.
    Both ``A`` and the elements of ``C`` have context object ``ctxAC``, while ``B`` has context object ``ctxB``.
    The length of the array ``C`` is the number of variables in ``ctxB``.
    Neither ``A`` nor ``B`` is allowed to alias any other polynomial.
    Return `1` for success and `0` for failure.
    The main method attemps to perform the calculation using matrices and chooses heuristically between the ``geobucket`` and ``horner`` methods if needed.

.. function:: void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A, const fmpz_mpoly_t B, const slong * c, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

    Set ``A`` to the evaluation of ``B`` where the variable of index ``i`` in ``ctxB`` is replaced by the variable of index ``c[i]`` in ``ctxAC``.
    The length of the array ``C`` is the number of variables in ``ctxB``.
    If any ``c[i]`` is negative, the corresponding variable of ``B`` is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in ``ctxAC``.


Multiplication
--------------------------------------------------------------------------------


.. function:: void fmpz_mpoly_mul(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C``.

.. function:: void fmpz_mpoly_mul_johnson(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` times ``C`` using Johnson's heap-based method.
    The first version always uses one thread.

.. function:: int fmpz_mpoly_mul_array(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
              int fmpz_mpoly_mul_array_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Try to set ``A`` to ``B`` times ``C`` using arrays.
    If the return is ``0``, the operation was unsuccessful. Otherwise, it was successful and the return is ``1``.
    The first version always uses one thread.

.. function:: int fmpz_mpoly_mul_dense(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    Try to set ``A`` to ``B`` times ``C`` using dense arithmetic.
    If the return is ``0``, the operation was unsuccessful. Otherwise, it was successful and the return is ``1``.


Powering
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fmpz_mpoly_pow_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t k, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` raised to the `k`-th power.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` raised to the `k`-th power.
    Return `1` for success and `0` for failure.


Division
--------------------------------------------------------------------------------

.. function:: int fmpz_mpoly_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    If ``A`` is divisible by ``B``, set ``Q`` to the exact quotient and return ``1``. Otherwise, set ``Q`` to zero and return ``0``.

.. function:: int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
              int fmpz_mpoly_divides_heap_threaded(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, slong thread_limit)

    Perform the operation of func::fmpz_mpoly_divides using the algorithm of Michael Monagan and Roman Pearce.
    The threaded version takes an upper limit on the number of threads to use, while the first version always uses one thread.

.. function:: void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Set ``Q`` and ``R`` to the quotient and remainder of ``A`` divided by ``B``. The monomials in ``R`` divisible by the leading monomial of ``B`` will have coefficients reduced modulo the absolute value of the leading coefficient of ``B``.
    Note that this function is not very useful if the leading coefficient ``B`` is not a unit.

.. function:: void fmpz_mpoly_quasidivrem(fmpz_t scale, fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Set ``scale``, ``Q`` and ``R`` so that ``Q`` and ``R` are the quotient and remainder of ``scale*A`` divided by ``B``. No monomials in ``R`` will be divisible by the leading monomial of ``B``.

.. function:: void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Perform the operation of :func:`fmpz_mpoly_divrem` and discard ``R``.
    Note that this function is not very useful if the division is not exact and the leading coefficient ``B`` is not a unit.

.. function:: void fmpz_mpoly_quasidiv(fmpz_t scale, fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Perform the operation of :func:`fmpz_mpoly_quasidivrem` and discard ``R``.

.. function:: void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx)

    This function is as per :func:`fmpz_mpoly_divrem` except that it takes an array of divisor polynomials ``B`` and it returns an array of quotient polynomials ``Q``.
    The number of divisor (and hence quotient) polynomials, is given by ``len``.
    Note that this function is not very useful if there is no unit among the leading coefficients in the array ``B``.

.. function:: void fmpz_mpoly_quasidivrem_ideal(fmpz_t scale, fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx)

    This function is as per :func:`fmpz_mpoly_quasidivrem` except that it takes an array of divisor polynomials ``B`` and it returns an array of quotient polynomials ``Q``.
    The number of divisor (and hence quotient) polynomials, is given by ``len``.


Greatest Common Divisor
--------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    Set ``M`` to the GCD of the terms of ``A``.
    If ``A`` is zero, ``M`` will be zero. Otherwise, ``M`` will be a monomial with positive coefficient.

.. function:: int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Try to set ``G`` to the GCD of ``A`` and ``B`` with positive leading coefficient. The GCD of zero and zero is defined to be zero.
    If the return is ``1`` the function was successful. Otherwise the return is  ``0`` and ``G`` is left untouched.

.. function:: int fmpz_mpoly_gcd_cofactors(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Do the operation of :func:`fmpz_mpoly_gcd` and also compute ``Abar = A/G`` and ``Bbar = B/G`` if successful.

.. function:: int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
              int fmpz_mpoly_gcd_brown_threaded(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Try to set ``G`` to the GCD of ``A`` and ``B`` using Brown's algorithm.
    The first version always uses one thread.

.. function:: int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    Try to set ``G`` to the GCD of ``A`` and ``B`` using Zippel's interpolation algorithm to interpolate coefficients from univariate images in the most significant variable.


Univariate Functions
--------------------------------------------------------------------------------

    An ``fmpz_mpoly_univar_t`` holds a univariate polynomial in some main variable
    with ``fmpz_mpoly_t`` coefficients in the remaining variables. These functions
    are useful when one wants to rewrite an element of `\mathbb{Z}[x_1, \dots, x_m]`
    as an element of `(\mathbb{Z}[x_1, \dots, x_{v-1}, x_{v+1}, \dots, x_m])[x_v]`
    and vise versa.

.. function:: void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

    Initialize `A`.

.. function:: void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

    Clear `A`.

.. function:: void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t A, fmpz_mpoly_univar_t B, const fmpz_mpoly_ctx_t ctx)

    Swap `A` and `B`.

.. function:: void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to a univariate form of ``B`` by pulling out the variable of index ``var``.
    The coefficients of ``A`` will still belong to the content ``ctx`` but will not depend on the variable of index ``var``.

.. function:: void fmpz_mpoly_from_univar(fmpz_mpoly_t A, const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to the normal form of ``B`` by putting in the variable of index ``var``.
    This function is undefined if the coefficients of ``B`` depend on the variable of index ``var``.

.. function:: int fmpz_mpoly_univar_degree_fits_si(const fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

    Return `1` if the degree of ``A`` with respect to the main variable fits an ``slong``. Otherwise, return `0`.

.. function:: slong fmpz_mpoly_univar_length(const fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

    Return the number of terms in ``A`` with respect to the main variable.

.. function:: slong fmpz_mpoly_univar_get_term_exp_si(fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Return the exponent of the term of index ``i`` of ``A``.

.. function:: void fmpz_mpoly_univar_get_term_coeff(fmpz_mpoly_t c, const fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_univar_swap_term_coeff(fmpz_mpoly_t c, fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    Set (resp. swap) ``c`` to (resp. with) the coefficient of the term of index ``i`` of ``A``.


Internal Functions
--------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_inflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx)

    Apply the function ``e -> shift[v] + stride[v]*e`` to each exponent ``e`` corresponding to the variable ``v``.
    It is assumed that each shift and stride is not negative.

.. function:: void fmpz_mpoly_deflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx)

    Apply the function ``e -> (e - shift[v])/stride[v]`` to each exponent ``e`` corresponding to the variable ``v``.
    If any ``stride[v]`` is zero, the corresponding numerator ``e - shift[v]`` is assumed to be zero, and the quotient is defined as zero.
    This allows the function to undo the operation performed by :func:`fmpz_mpoly_inflate` when possible.

.. function:: void fmpz_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    For each variable `v` let `S_v` be the set of exponents appearing on `v`.
    Set ``shift[v]`` to `\operatorname{min}(S_v)` and set ``stride[v]`` to `\operatorname{gcd}(S-\operatorname{min}(S_v))`.
    If ``A`` is zero, all shifts and strides are set to zero.


.. function:: void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx)

    Set ``A`` to ``B`` raised to the `k`-th power, using the Monagan and Pearce FPS algorithm.
    It is assumed that ``B`` is not zero and `k \geq 2`.

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

