.. _fmpz-mod-mpoly:

**fmpz_mod_mpoly.h** -- polynomials over the integers mod n
===============================================================================

    The exponents follow the ``mpoly`` interface.
    A coefficient may be referenced as a ``fmpz *``, but this may disappear in a future version.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mod_mpoly_struct

    A structure holding a multivariate polynomial over the integers mod n.

.. type:: fmpz_mod_mpoly_t

    An array of length `1` of ``fmpz_mod_mpoly_ctx_struct``.

.. type:: fmpz_mod_mpoly_ctx_struct

    Context structure representing the parent ring of an ``fmpz_mod_mpoly``.

.. type:: fmpz_mod_mpoly_ctx_t

    An array of length `1` of ``fmpz_mod_mpoly_struct``.


Context object
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz_t p)

    Initialise a context object for a polynomial ring modulo *n* with *nvars* variables and ordering *ord*.
    The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

.. function:: slong fmpz_mod_mpoly_ctx_nvars(fmpz_mod_mpoly_ctx_t ctx)

    Return the number of variables used to initialize the context.

.. function:: ordering_t fmpz_mod_mpoly_ctx_ord(const fmpz_mod_mpoly_ctx_t ctx)

    Return the ordering used to initialize the context.

.. function:: void fmpz_mod_mpoly_ctx_get_modulus(fmpz_t n, const fmpz_mod_mpoly_ctx_t ctx)

    Set *n* to the modulus used to initialize the context.

.. function:: void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx)

    Release up any space allocated by an *ctx*.


Memory management
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_init(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Initialise *A* for use with the given an initialised context object. Its value is set to zero.

.. function:: void fmpz_mod_mpoly_init2(fmpz_mod_mpoly_t A, slong alloc, const fmpz_mod_mpoly_ctx_t ctx)

    Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for *alloc* terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

.. function:: void fmpz_mod_mpoly_init3(fmpz_mod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)

    Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    It is allocated with space for *alloc* terms and *bits* bits for the exponents.

.. function:: void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Release any space allocated for *A*.


Input/Output
--------------------------------------------------------------------------------

    The variable strings in *x* start with the variable of most significance at index `0`. If *x* is ``NULL``, the variables are named ``x1``, ``x2``, ect.

.. function:: char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing *A*, given an array of variable strings *x*.

.. function:: int fmpz_mod_mpoly_fprint_pretty(FILE * file, const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Print a string representing *A* to *file*.

.. function:: int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Print a string representing *A* to ``stdout``.

.. function:: int fmpz_mod_mpoly_set_str_pretty(fmpz_mod_mpoly_t A, const char * str, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the polynomial in the null-terminates string *str* given an array *x* of variable strings.
    If parsing *str* fails, *A* is set to zero, and `-1` is returned. Otherwise, `0`  is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*. The character ``^`` must be immediately followed by the (integer) exponent.
    If any division is not exact, parsing fails.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the variable of index *var*, where `var = 0` corresponds to the variable with the most significance with respect to the ordering. 

.. function:: int fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    If `var \ge 0`, return `1` if *A* is equal to the `var`-th generator, otherwise return `0`.
    If `var < 0`, return `1` if the polynomial is equal to any generator, otherwise return `0`.

.. function:: void fmpz_mod_mpoly_set(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    
    Set *A* to *B*.

.. function:: int fmpz_mod_mpoly_equal(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is equal to *B*, else return `0`.

.. function:: void fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t poly1, fmpz_mod_mpoly_t poly2, const fmpz_mod_mpoly_ctx_t ctx)

    Efficiently swap *A* and *B*.


Constants
--------------------------------------------------------------------------------


.. function:: int fmpz_mod_mpoly_is_fmpz(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is a constant, else return `0`.

.. function:: void fmpz_mod_mpoly_get_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Assuming that *A* is a constant, set *c* to this constant.
    This function throws if *A* is not a constant.

.. function:: void fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_ui(fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the constant *c*.

.. function:: void fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the constant `0`.

.. function:: void fmpz_mod_mpoly_one(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the constant `1`.

.. function:: int fmpz_mod_mpoly_equal_fmpz(const fmpz_mod_mpoly_t A, fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_equal_ui(const fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_equal_si(const fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is equal to the constant *c*, else return `0`.

.. function:: int fmpz_mod_mpoly_is_zero(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is the constant `0`, else return `0`.

.. function:: int fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is the constant `1`, else return `0`.


Degrees
--------------------------------------------------------------------------------


.. function:: int fmpz_mod_mpoly_degrees_fit_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if the degrees of *A* with respect to each variable fit into an ``slong``, otherwise return `0`.

.. function:: void fmpz_mod_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_degrees_si(slong * degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *degs* to the degrees of *A* with respect to each variable.
    If *A* is zero, all degrees are set to `-1`.

.. function:: void fmpz_mod_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
              slong fmpz_mod_mpoly_degree_si(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Either return or set *deg* to the degree of *A* with respect to the variable of index *var*.
    If *A* is zero, the degree is defined to be `-1`.

.. function:: int fmpz_mod_mpoly_total_degree_fits_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if the total degree of *A* fits into an ``slong``, otherwise return `0`.

.. function:: void fmpz_mod_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
              slong fmpz_mod_mpoly_total_degree_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Either return or set *tdeg* to the total degree of *A*.
    If *A* is zero, the total degree is defined to be `-1`.

.. function:: void fmpz_mod_mpoly_used_vars(int * used, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    For each variable index *i*, set ``used[i]`` to nonzero if the variable of index *i* appears in *A* and to zero otherwise.


Coefficients
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)

    Assuming that *M* is a monomial, set *c* to the coefficient of the corresponding monomial in *A*.
    This function thows if *M* is not a monomial.

.. function:: void fmpz_mod_mpoly_set_coeff_fmpz_monomial(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)

    Assuming that *M* is a monomial, set the coefficient of the corresponding monomial in *A* to *c*.
    This function thows if *M* is not a monomial.

.. function:: void fmpz_mod_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mod_mpoly_t A, ulong const * exp, const fmpz_mod_mpoly_ctx_t ctx)

    Set *c* to the coefficient of the monomial with exponent vector *exp*.

.. function:: void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_coeff_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_coeff_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, ulong const * exp, fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_coeff_ui_ui(fmpz_mod_mpoly_t A, ulong c, ulong const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_coeff_si_ui(fmpz_mod_mpoly_t A, slong c, ulong const * exp, const fmpz_mod_mpoly_ctx_t ctx)

    Set the coefficient of the monomial with exponent vector *exp* to *c*.

.. function:: void fmpz_mod_mpoly_get_coeff_vars_ui(fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpz_mod_mpoly_ctx_t ctx)

    Set *C* to the coefficient of *A* with respect to the variables in *vars* with powers in the corresponding array *exps*.
    Both *vars* and *exps* point to array of length *length*. It is assumed that `0 < length \le nvars(A)` and that the variables in *vars* are distinct.


Comparison
--------------------------------------------------------------------------------


.. function:: int fmpz_mod_mpoly_cmp(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` (resp. `-1`, or `0`) if *A* is after (resp. before, same as) *B* in some arbitrary but fixed total ordering of the polynomials.
    This ordering agrees with the usual ordering of monomials when *A* and *B* are both monomials.


Container operations
--------------------------------------------------------------------------------

    These functions deal with violations of the internal canonical representation.
    If a term index is negative or not strictly less than the length of the polynomial, the function will throw.

.. function:: int fmpz_mod_mpoly_is_canonical(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is in canonical form. Otherwise, return `0`.
    To be in canonical form, all of the terms must have nonzero coefficient, and the terms must be sorted from greatest to least.

.. function:: slong fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return the number of terms in *A*.
    If the polynomial is in canonical form, this will be the number of nonzero coefficients.

.. function:: void fmpz_mod_mpoly_resize(fmpz_mod_mpoly_t A, slong new_length, const fmpz_mod_mpoly_ctx_t ctx)

    Set the length of *A* to ``new_length``.
    Terms are either deleted from the end, or new zero terms are appended.

.. function:: void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set *c* to the coefficient of the term of index *i*.

.. function:: void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A, slong i, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_term_coeff_ui(fmpz_mod_mpoly_t A, slong i, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_term_coeff_si(fmpz_mod_mpoly_t A, slong i, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Set the coefficient of the term of index *i* to *c*.

.. function:: int fmpz_mod_mpoly_term_exp_fits_si(const fmpz_mod_mpoly_t poly, slong i, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_term_exp_fits_ui(const fmpz_mod_mpoly_t poly, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if all entries of the exponent vector of the term of index *i* fit into an ``slong`` (resp. a ``ulong``). Otherwise, return `0`.

.. function:: void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_get_term_exp_si(slong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set *exp* to the exponent vector of the term of index *i*.
    The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

.. function:: ulong fmpz_mod_mpoly_get_term_var_exp_ui(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)
              slong fmpz_mod_mpoly_get_term_var_exp_si(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Return the exponent of the variable *var* of the term of index *i*.
    This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

.. function:: void fmpz_mod_mpoly_set_term_exp_fmpz(fmpz_mod_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, slong i, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)

    Set the exponent vector of the term of index *i* to *exp*.

.. function:: void fmpz_mod_mpoly_get_term(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set *M* to the term of index *i* in *A*.

.. function:: void fmpz_mod_mpoly_get_term_monomial(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set *M* to the monomial of the term of index *i* in *A*. The coefficient of *M* will be one.

.. function:: void fmpz_mod_mpoly_push_term_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_push_term_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_push_term_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_push_term_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_push_term_ui_ui(fmpz_mod_mpoly_t A, ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_push_term_si_ui(fmpz_mod_mpoly_t A, slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)

    Append a term to *A* with coefficient *c* and exponent vector *exp*.
    This function runs in constant average time.

.. function:: void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Sort the terms of *A* into the canonical ordering dictated by the ordering in *ctx*.
    This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    This function runs in linear time in the size of *A*.

.. function:: void fmpz_mod_mpoly_combine_like_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Combine adjacent like terms in *A* and delete terms with coefficient zero.
    If the terms of *A* were sorted to begin with, the result will be in canonical form.
    This function runs in linear time in the size of *A*.

.. function:: void fmpz_mod_mpoly_reverse(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the reversal of *B*.


Random generation
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bound - 1]``.
    The exponents of each variable are generated by calls to ``n_randint(state, exp_bound)``.

.. function:: void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const fmpz_mod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bounds[i] - 1]``.
    The exponents of the variable of index *i* are generated by calls to ``n_randint(state, exp_bounds[i])``.

.. function:: void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const fmpz_mod_mpoly_ctx_t ctx)

    Generate a random polynomial with length up to *length* and exponents whose packed form does not exceed the given bit count.


Addition/Subtraction
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_add_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_add_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_add_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B + c`.

.. function:: void fmpz_mod_mpoly_sub_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_sub_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_sub_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B - c`.

.. function:: void fmpz_mod_mpoly_add(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B + C`.

.. function:: void fmpz_mod_mpoly_sub(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B - C`.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_neg(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `-B`.

.. function:: void fmpz_mod_mpoly_scalar_mul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_scalar_mul_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_scalar_mul_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B \times c`.

.. function:: void fmpz_mod_mpoly_scalar_addmul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_t d, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *A* to `B + C \times d`.

.. function:: void fmpz_mod_mpoly_make_monic(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to *B* divided by the leading coefficient of *B*. This throws if *B* is zero or the leading coefficient is not invertible.


Differentiation
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_derivative(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the derivative of *B* with respect to the variable of index *var*.


Evaluation
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: void fmpz_mod_mpoly_evaluate_all_fmpz(fmpz_t eval, const fmpz_mod_mpoly_t A, fmpz * const * vals, const fmpz_mod_mpoly_ctx_t ctx)

    Set *ev* to the evaluation of *A* where the variables are replaced by the corresponding elements of the array *vals*.

.. function:: void fmpz_mod_mpoly_evaluate_one_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_t val, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the evaluation of *B* where the variable of index *var* is replaced by *val*.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mod_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mod_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB)

    Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    The context object of *B* is *ctxB*.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_geobucket(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
              int fmpz_mod_mpoly_compose_fmpz_mod_mpoly(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)

    Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    Both *A* and the elements of *C* have context object *ctxAC*, while *B* has context object *ctxB*.
    The length of the array *C* is the number of variables in *ctxB*.
    Neither *A* nor *B* is allowed to alias any other polynomial.
    Return `1` for success and `0` for failure.
    The main method attempts to perform the calculation using matrices and chooses heuristically between the ``geobucket`` and ``horner`` methods if needed.

.. function:: void fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const slong * c, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)

    Set *A* to the evaluation of *B* where the variable of index *i* in *ctxB* is replaced by the variable of index ``c[i]`` in *ctxAC*.
    The length of the array *C* is the number of variables in *ctxB*.
    If any ``c[i]`` is negative, the corresponding variable of *B* is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in *ctxAC*.


Multiplication
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_mpoly_mul(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B \times C`.

.. function:: void fmpz_mod_mpoly_mul_johnson(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to `B \times C` using Johnson's heap-based method.

.. function:: int fmpz_mod_mpoly_mul_dense(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *A* to `B \times C` using dense arithmetic.
    If the return is `0`, the operation was unsuccessful. Otherwise, it was successful and the return is `1`.


Powering
--------------------------------------------------------------------------------

    These functions return `0` when the operation would imply unreasonable arithmetic.

.. function:: int fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t k, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to *B* raised to the `k`-th power.
    Return `1` for success and `0` for failure.

.. function:: int fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to *B* raised to the `k`-th power.
    Return `1` for success and `0` for failure.


Division
--------------------------------------------------------------------------------

The division functions assume that the modulus is prime.

.. function:: int fmpz_mod_mpoly_divides(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    If *A* is divisible by *B*, set *Q* to the exact quotient and return `1`. Otherwise, set *Q* to zero and return `0`.

.. function:: void fmpz_mod_mpoly_div(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Set *Q* to the quotient of *A* by *B*, discarding the remainder.

.. function:: void fmpz_mod_mpoly_divrem(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Set *Q* and *R* to the quotient and remainder of *A* divided by *B*.

.. function:: void fmpz_mod_mpoly_divrem_ideal(fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len, const fmpz_mod_mpoly_ctx_t ctx)

    This function is as per :func:`fmpz_mod_mpoly_divrem` except that it takes an array of divisor polynomials *B* and it returns an array of quotient polynomials *Q*.
    The number of divisor (and hence quotient) polynomials, is given by *len*.


Greatest Common Divisor
--------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_term_content(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Set *M* to the GCD of the terms of *A*.
    If *A* is zero, *M* will be zero. Otherwise, *M* will be a monomial with coefficient one.

.. function:: int fmpz_mod_mpoly_content_vars(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_t A, slong * vars, slong vars_length, const fmpz_mod_mpoly_ctx_t ctx)

    Set *g* to the GCD of the cofficients of *A* when viewed as a polynomial in the variables *vars*.
    Return `1` for success and `0` for failure. Upon succcess, *g* will be independent of the variables *vars*.

.. function:: int fmpz_mod_mpoly_gcd(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *G* to the monic GCD of *A* and *B*. The GCD of zero and zero is defined to be zero.
    If the return is `1` the function was successful. Otherwise the return is `0` and *G* is left untouched.

.. function:: int fmpz_mod_mpoly_gcd_cofactors(fmpz_mod_mpoly_t G, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Do the operation of :func:`fmpz_mod_mpoly_gcd` and also compute `Abar = A/G` and `Bbar = B/G` if successful.

.. function:: int fmpz_mod_mpoly_gcd_brown(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_gcd_hensel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_gcd_subresultant(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_gcd_zippel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_gcd_zippel2(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *G* to the GCD of *A* and *B* using various algorithms.

.. function:: int fmpz_mod_mpoly_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *R* to the resultant of *A* and *B* with respect to the variable of index *var*.

.. function:: int fmpz_mod_mpoly_discriminant(fmpz_mod_mpoly_t D, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *D* to the discriminant of *A* with respect to the variable of index *var*.


Square Root
--------------------------------------------------------------------------------

The square root functions assume that the modulus is prime for correct operation.

.. function:: int fmpz_mod_mpoly_sqrt(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    If `Q^2=A` has a solution, set *Q* to a solution and return `1`, otherwise return `0` and set *Q* to zero.

.. function:: int fmpz_mod_mpoly_is_square(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if *A* is a perfect square, otherwise return `0`.

.. function:: int fmpz_mod_mpoly_quadratic_root(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

    If `Q^2+AQ=B` has a solution, set *Q* to a solution and return `1`, otherwise return `0`.


Univariate Functions
--------------------------------------------------------------------------------

    An ``fmpz_mod_mpoly_univar_t`` holds a univariate polynomial in some main variable
    with ``fmpz_mod_mpoly_t`` coefficients in the remaining variables. These functions
    are useful when one wants to rewrite an element of `\mathbb{Z}/n\mathbb{Z}[x_1, \dots, x_m]`
    as an element of `(\mathbb{Z}/n\mathbb{Z}[x_1, \dots, x_{v-1}, x_{v+1}, \dots, x_m])[x_v]`
    and vise versa.

.. function:: void fmpz_mod_mpoly_univar_init(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Initialize *A*.

.. function:: void fmpz_mod_mpoly_univar_clear(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Clear *A*.

.. function:: void fmpz_mod_mpoly_univar_swap(fmpz_mod_mpoly_univar_t A, fmpz_mod_mpoly_univar_t B, const fmpz_mod_mpoly_ctx_t ctx)

    Swap *A* and *B*.

.. function:: void fmpz_mod_mpoly_to_univar(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to a univariate form of *B* by pulling out the variable of index *var*.
    The coefficients of *A* will still belong to the content *ctx* but will not depend on the variable of index *var*.

.. function:: void fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_univar_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

    Set *A* to the normal form of *B* by putting in the variable of index *var*.
    This function is undefined if the coefficients of *B* depend on the variable of index *var*.

.. function:: int fmpz_mod_mpoly_univar_degree_fits_si(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return `1` if the degree of *A* with respect to the main variable fits an ``slong``. Otherwise, return `0`.

.. function:: slong fmpz_mod_mpoly_univar_length(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)

    Return the number of terms in *A* with respect to the main variable.

.. function:: slong fmpz_mod_mpoly_univar_get_term_exp_si(fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Return the exponent of the term of index *i* of *A*.

.. function:: void fmpz_mod_mpoly_univar_get_term_coeff(fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_univar_swap_term_coeff(fmpz_mod_mpoly_t c, fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Set (resp. swap) *c* to (resp. with) the coefficient of the term of index *i* of *A*.

.. function:: void fmpz_mod_mpoly_univar_set_coeff_ui(fmpz_mod_mpoly_univar_t Ax, ulong e, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)

    Set the coefficient of `X^e` in *Ax* to *c*.

.. function:: int fmpz_mod_mpoly_univar_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *R* to the resultant of *Ax* and *Bx*.

.. function:: int fmpz_mod_mpoly_univar_discriminant(fmpz_mod_mpoly_t D, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_ctx_t ctx)

    Try to set *D* to the discriminant of *Ax*.


Internal Functions
--------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_inflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)

    Apply the function ``e -> shift[v] + stride[v]*e`` to each exponent ``e`` corresponding to the variable ``v``.
    It is assumed that each shift and stride is not negative.

.. function:: void fmpz_mod_mpoly_deflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)

    Apply the function ``e -> (e - shift[v])/stride[v]`` to each exponent ``e`` corresponding to the variable ``v``.
    If any ``stride[v]`` is zero, the corresponding numerator ``e - shift[v]`` is assumed to be zero, and the quotient is defined as zero.
    This allows the function to undo the operation performed by :func:`fmpz_mod_mpoly_inflate` when possible.

.. function:: void fmpz_mod_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

    For each variable `v` let `S_v` be the set of exponents appearing on `v`.
    Set ``shift[v]`` to `\operatorname{min}(S_v)` and set ``stride[v]`` to `\operatorname{gcd}(S-\operatorname{min}(S_v))`.
    If *A* is zero, all shifts and strides are set to zero.

