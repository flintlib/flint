.. _fmpz-mod-poly-factor:

**fmpz_mod_poly_factor.h** -- factorisation of polynomials over integers mod n
==================================================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mod_poly_factor_struct

    A structure representing a polynomial in factorised form as a product of
    polynomials with associated exponents.

.. type:: fmpz_mod_poly_factor_t

    An array of length 1 of ``fmpz_mpoly_factor_struct``.

Factorisation
--------------------------------------------------------------------------------

The factorization, irreducibility testing and root finding functions in
this section are wrappers around the generic implementations in the
``gr_poly`` module (see :func:`gr_poly_factor_finite_field`,
:func:`gr_poly_is_irreducible` and :func:`gr_poly_roots_finite_field`),
which select algorithms and cutoffs internally and use several threads
when these are available.


The factorisation, irreducibility and squarefreeness functions in this
module assume that the modulus is prime. This is not checked.


.. function:: void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

    Initialises ``fac`` for use.

.. function:: void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

    Frees all memory associated with ``fac``.

.. function:: void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac, slong alloc, const fmpz_mod_ctx_t ctx)

    Reallocates the factor structure to provide space for
    precisely ``alloc`` factors.

.. function:: void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac, slong len, const fmpz_mod_ctx_t ctx)

    Ensures that the factor structure has space for at
    least ``len`` factors.  This function takes care
    of the case of repeated calls by always at least
    doubling the number of factors the structure can hold.

.. function:: void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

    Sets ``res`` to the same factorisation as ``fac``.

.. function:: void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

    Prints the entries of ``fac`` to standard output.

.. function:: void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac, const fmpz_mod_poly_t poly, slong exp, const fmpz_mod_ctx_t ctx)

    Inserts the factor ``poly`` with multiplicity ``exp`` into
    the factorisation ``fac``.

    If ``fac`` already contains ``poly``, then ``exp`` simply
    gets added to the exponent of the existing entry.

.. function:: void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac, const fmpz_mod_ctx_t ctx)

    Concatenates two factorisations.

    This is equivalent to calling :func:`fmpz_mod_poly_factor_insert`
    repeatedly with the individual factors of ``fac``.

    Does not support aliasing between ``res`` and ``fac``.

.. function:: void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp, const fmpz_mod_ctx_t ctx)

    Raises ``fac`` to the power ``exp``.

.. function:: int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.

.. function:: int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses fast distinct-degree factorisation.

.. function:: int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses Rabin irreducibility test.

.. function:: int fmpz_mod_poly_is_irreducible_rabin_f(fmpz_t r, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    Either sets `r` to `1` and returns 1 if the polynomial ``f`` is
    irreducible or `0` otherwise, or sets `r` to a nontrivial factor of
    `p`.

    This algorithm correctly determines whether `f` is irreducible over
    `\mathbb{Z}/p\mathbb{Z}`, even for composite `f`, or it finds a factor
    of `p`.

.. function:: int _fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_mod_ctx_t ctx)

    Returns 1 if ``(f, len)`` is squarefree, and 0 otherwise. As a
    special case, the zero polynomial is not considered squarefree.
    There are no restrictions on the length.

.. function:: int _fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz * f, slong len, const fmpz_mod_ctx_t ctx)

    If `fac` returns with the value `1` then the function operates as per
    :func:`_fmpz_mod_poly_is_squarefree`, otherwise `f` is set to a nontrivial
    factor of `p`.

.. function:: int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    Returns 1 if ``f`` is squarefree, and 0 otherwise. As a special
    case, the zero polynomial is not considered squarefree.

.. function:: int fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    If `fac` returns with the value `1` then the function operates as per
    :func:`fmpz_mod_poly_is_squarefree`, otherwise `f` is set to a nontrivial
    factor of `p`.

Root Finding
--------------------------------------------------------------------------------

.. function:: void fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f, int with_multiplicity, const fmpz_mod_ctx_t ctx)

    Fill `r` with factors of the form `x - r_i` where the `r_i` are the distinct roots of a nonzero `f` in `Z/pZ`.
    It is expected and not checked that the modulus of `ctx` is prime.
    If `with\_multiplicity` is zero, the exponent `e_i` of the factor `x - r_i` is `1`. Otherwise, it is the largest `e_i` such that `(x-r_i)^e_i` divides `f`.
    This function throws if `f` is zero, but is otherwise always successful.

.. function:: int fmpz_mod_poly_roots_factored(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f, int with_multiplicity, const fmpz_factor_t n, const fmpz_mod_ctx_t ctx)

    Fill `r` with factors of the form `x - r_i` where the `r_i` are the distinct roots of a nonzero `f` in `Z/nZ`.
    It is expected and not checked that `n` is a prime factorization of the modulus of `ctx`.
    If `with\_multiplicity` is zero, the exponent `e_i` of the factor `x - r_i` is `1`. Otherwise, it is the largest `e_i` such that `(x-r_i)^e_i` divides `f`.
    The roots are first found modulo the primes in `n`, then lifted to the corresponding prime powers, then combined into roots of the original polynomial `f`.
    A return of `1` indicates the function was successful. A return of `0` indicates the function was not able to find the roots, possibly because there are too many of them.
    This function throws if `f` is zero.
