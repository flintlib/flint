.. _fq-poly-factor:

**fq_poly_factor.h** -- factorisation of univariate polynomials over finite fields
==================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_poly_factor_struct

.. type:: fq_poly_factor_t

Memory Management
--------------------------------------------------------------------------------


.. function:: void fq_poly_factor_init(fq_poly_factor_t fac, const fq_ctx_t ctx)

    Initialises ``fac`` for use. An :type:`fq_poly_factor_t`
    represents a polynomial in factorised form as a product of
    polynomials with associated exponents.

.. function:: void fq_poly_factor_clear(fq_poly_factor_t fac, const fq_ctx_t ctx)

    Frees all memory associated with ``fac``.

.. function:: void fq_poly_factor_realloc(fq_poly_factor_t fac, slong alloc, const fq_ctx_t ctx)

    Reallocates the factor structure to provide space for
    precisely ``alloc`` factors.

.. function:: void fq_poly_factor_fit_length(fq_poly_factor_t fac, slong len, const fq_ctx_t ctx)

    Ensures that the factor structure has space for at least
    ``len`` factors.  This function takes care of the case of
    repeated calls by always at least doubling the number of factors
    the structure can hold.


Basic Operations
--------------------------------------------------------------------------------


.. function:: void fq_poly_factor_set(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)

    Sets ``res`` to the same factorisation as ``fac``.

.. function:: void fq_poly_factor_print_pretty(const fq_poly_factor_t fac, const char * var, const fq_ctx_t ctx)

    Pretty-prints the entries of ``fac`` to standard output.

.. function:: void fq_poly_factor_print(const fq_poly_factor_t fac, const fq_ctx_t ctx)

    Prints the entries of ``fac`` to standard output.

.. function:: void fq_poly_factor_insert(fq_poly_factor_t fac, const fq_poly_t poly, slong exp, const fq_ctx_t ctx)

    Inserts the factor ``poly`` with multiplicity ``exp`` into
    the factorisation ``fac``.

    If ``fac`` already contains ``poly``, then ``exp`` simply
    gets added to the exponent of the existing entry.

.. function:: void fq_poly_factor_concat(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)

    Concatenates two factorisations.

    This is equivalent to calling :func:`fq_poly_factor_insert`
    repeatedly with the individual factors of ``fac``.

    Does not support aliasing between ``res`` and ``fac``.

.. function:: void fq_poly_factor_pow(fq_poly_factor_t fac, slong exp, const fq_ctx_t ctx)

    Raises ``fac`` to the power ``exp``.

.. function:: ulong fq_poly_remove(fq_poly_t f, const fq_poly_t p, const fq_ctx_t ctx)

    Removes the highest possible power of ``p`` from ``f`` and
    returns the exponent.


Irreducibility Testing
--------------------------------------------------------------------------------

.. function:: int fq_poly_is_irreducible(const fq_poly_t f, const fq_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.

.. function:: int fq_poly_is_irreducible_ddf(const fq_poly_t f, const fq_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses fast distinct-degree factorisation.

.. function:: int fq_poly_is_irreducible_ben_or(const fq_poly_t f, const fq_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses Ben-Or's irreducibility test.

.. function:: int _fq_poly_is_squarefree(const fq_struct * f, slong len, const fq_ctx_t ctx)

    Returns 1 if ``(f, len)`` is squarefree, and 0 otherwise. As a
    special case, the zero polynomial is not considered squarefree.
    There are no restrictions on the length.

.. function:: int fq_poly_is_squarefree(const fq_poly_t f, const fq_ctx_t ctx)

    Returns 1 if ``f`` is squarefree, and 0 otherwise. As a special
    case, the zero polynomial is not considered squarefree.



Factorisation
--------------------------------------------------------------------------------

The factorization, irreducibility testing and root finding functions in
this section are wrappers around the generic implementations in the
``gr_poly`` module (see :func:`gr_poly_factor_finite_field`,
:func:`gr_poly_is_irreducible` and :func:`gr_poly_roots_finite_field`),
which select algorithms and cutoffs internally and use several threads
when these are available.



Root Finding
--------------------------------------------------------------------------------

.. function:: void fq_poly_roots(fq_poly_factor_t r, const fq_poly_t f, int with_multiplicity, const fq_ctx_t ctx)

    Fill `r` with factors of the form `x - r_i` where the `r_i` are the distinct roots of a nonzero `f` in `F_q`.
    If `with\_multiplicity` is zero, the exponent `e_i` of the factor `x - r_i` is `1`. Otherwise, it is the largest `e_i` such that `(x-r_i)^e_i` divides `f`.
    This function throws if `f` is zero, but is otherwise always successful.
