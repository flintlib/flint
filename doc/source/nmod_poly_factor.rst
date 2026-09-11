.. _nmod-poly-factor:

**nmod_poly_factor.h** -- factorisation of univariate polynomials over integers mod n (word-size n)
===================================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: nmod_poly_factor_struct

.. type:: nmod_poly_factor_t

Factorisation
--------------------------------------------------------------------------------

The factorization, irreducibility testing and root finding functions in
this module are wrappers around the generic implementations in
the ``gr_poly`` module (see :func:`gr_poly_factor_finite_field`,
:func:`gr_poly_is_irreducible` and :func:`gr_poly_roots_finite_field`),
which select algorithms and cutoffs internally and use several threads
when available.


The factorisation and irreducibility functions in this module assume
that the modulus is prime. This is not checked.

.. function:: void nmod_poly_factor_init(nmod_poly_factor_t fac)

    Initialises ``fac`` for use. An ``nmod_poly_factor_t``
    represents a polynomial in factorised form as a product of
    polynomials with associated exponents.

.. function:: void nmod_poly_factor_clear(nmod_poly_factor_t fac)

    Frees all memory associated with ``fac``.

.. function:: void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc)

    Reallocates the factor structure to provide space for
    precisely ``alloc`` factors.

.. function:: void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, slong len)

    Ensures that the factor structure has space for at
    least ``len`` factors.  This function takes care
    of the case of repeated calls by always at least
    doubling the number of factors the structure can hold.

.. function:: void nmod_poly_factor_set(nmod_poly_factor_t res, const nmod_poly_factor_t fac)

    Sets ``res`` to the same factorisation as ``fac``.

.. function:: void nmod_poly_factor_print(const nmod_poly_factor_t fac)

    Prints the entries of ``fac`` to standard output.

.. function:: void nmod_poly_factor_insert(nmod_poly_factor_t fac, const nmod_poly_t poly, slong exp)

    Inserts the factor ``poly`` with multiplicity ``exp`` into
    the factorisation ``fac``.

    If ``fac`` already contains ``poly``, then ``exp`` simply
    gets added to the exponent of the existing entry.

.. function:: void nmod_poly_factor_concat(nmod_poly_factor_t res, const nmod_poly_factor_t fac)

    Concatenates two factorisations.

    This is equivalent to calling :func:`nmod_poly_factor_insert`
    repeatedly with the individual factors of ``fac``.

    Does not support aliasing between ``res`` and ``fac``.

.. function:: void nmod_poly_factor_pow(nmod_poly_factor_t fac, slong exp)

    Raises ``fac`` to the power ``exp``.

.. function:: int nmod_poly_is_irreducible(const nmod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.

.. function:: int nmod_poly_is_irreducible_ddf(const nmod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses fast distinct-degree factorisation.

.. function:: int nmod_poly_is_irreducible_rabin(const nmod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses Rabin irreducibility test.

.. function:: int _nmod_poly_is_squarefree(nn_srcptr f, slong len, nmod_t mod)

    Returns 1 if ``(f, len)`` is squarefree, and 0 otherwise. As a
    special case, the zero polynomial is not considered squarefree.
    There are no restrictions on the length.

.. function:: int nmod_poly_is_squarefree(const nmod_poly_t f)

    Returns 1 if ``f`` is squarefree, and 0 otherwise. As a special
    case, the zero polynomial is not considered squarefree.

.. function:: void nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)

    Sets ``res`` to a square-free factorization of ``f``.
