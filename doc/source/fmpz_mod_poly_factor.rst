.. _fmpz-mod-poly-factor:

**fmpz_mod_poly_factor.h** -- factorisation of polynomials over integers mod n
==================================================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mod_poly_factor_struct

.. type:: fmpz_mod_poly_factor_t

    Description.

Factorisation
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac)

    Initialises ``fac`` for use. An ``fmpz_mod_poly_factor_t``
    represents a polynomial in factorised form as a product of
    polynomials with associated exponents.

.. function:: void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac)

    Frees all memory associated with ``fac``.

.. function:: void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac, slong alloc)

    Reallocates the factor structure to provide space for
    precisely ``alloc`` factors.

.. function:: void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac, slong len)

    Ensures that the factor structure has space for at
    least ``len`` factors.  This function takes care
    of the case of repeated calls by always, at least
    doubling the number of factors the structure can hold.

.. function:: void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac)

    Sets ``res`` to the same factorisation as ``fac``.

.. function:: void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac)

    Prints the entries of ``fac`` to standard output.

.. function:: void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac, const fmpz_mod_poly_t poly, slong exp)

    Inserts the factor ``poly`` with multiplicity ``exp`` into
    the factorisation ``fac``.

    If ``fac`` already contains ``poly``, then ``exp`` simply
    gets added to the exponent of the existing entry.

.. function:: void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac)

    Concatenates two factorisations.

    This is equivalent to calling ``fmpz_mod_poly_factor_insert()``
    repeatedly with the individual factors of ``fac``.

    Does not support aliasing between ``res`` and ``fac``.

.. function:: void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp)

    Raises ``fac`` to the power ``exp``.

.. function:: int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.

.. function:: int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses fast distinct-degree factorisation.

.. function:: int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f)

    Returns 1 if the polynomial ``f`` is irreducible, otherwise returns 0.
    Uses Rabin irreducibility test.

.. function:: int fmpz_mod_poly_is_irreducible_rabin_f(fmpz_t f, const fmpz_mod_poly_t f)

    Either sets `f` to `1` and return 1 if the polynomial ``f`` is 
    irreducible or `0` otherwise, or set `f` to a nontrivial factor of
    `p`.

    This algorithm correctly determines whether `f` to is irreducible over 
    `\mathbb{Z}/p\mathbb{Z}`, even for composite `f`, or it finds a factor
    of `p`.

.. function:: int _fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_t p)

    Returns 1 if ``(f, len)`` is squarefree, and 0 otherwise. As a
    special case, the zero polynomial is not considered squarefree.
    There are no restrictions on the length.

.. function:: int _fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz * f, slong len, const fmpz_t p)

    If `fac` returns with the value `1` then the function operates as per
    ``_fmpz_mod_poly_is_squarefree``, otherwise `f` is set to a nontrivial
    factor of `p`.

.. function:: int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f)

    Returns 1 if ``f`` is squarefree, and 0 otherwise. As a special
    case, the zero polynomial is not considered squarefree.

.. function:: int fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz_mod_poly_t f)

    If `fac` returns with the value `1` then the function operates as per
    ``fmpz_mod_poly_is_squarefree``, otherwise `f` is set to a nontrivial
    factor of `p`.

.. function:: int fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor, flint_rand_t state, const fmpz_mod_poly_t pol, slong d)

    Probabilistic equal degree factorisation of ``pol`` into
    irreducible factors of degree ``d``. If it passes, a factor is
    placed in ``factor`` and 1 is returned, otherwise 0 is returned and
    the value of factor is undetermined.

    Requires that ``pol`` be monic, non-constant and squarefree.

.. function:: void fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors, const fmpz_mod_poly_t pol, slong d)

    Assuming ``pol`` is a product of irreducible factors all of
    degree ``d``, finds all those factors and places them in factors.
    Requires that ``pol`` be monic, non-constant and squarefree.

.. function:: void fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly, slong * const *degs)

    Factorises a monic non-constant squarefree polynomial ``poly``
    of degree n into factors `f[d]` such that for `1 \leq d \leq n`
    `f[d]` is the product of the monic irreducible factors of ``poly``
    of degree `d`. Factors `f[d]` are stored in ``res``, and the degree `d`
    of the irreducible factors is stored in ``degs`` in the same order
    as the factors.

    Requires that ``degs`` has enough space for `(n/2)+1 * sizeof(slong)`.

.. function:: void fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly, slong * const *degs)

    Multithreaded version of ``fmpz_mod_poly_factor_distinct_deg``.

.. function:: void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f)

    Sets ``res`` to a squarefree factorization of ``f``.

.. function:: void fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f)

    Factorises a non-constant polynomial ``f`` into monic irreducible
    factors choosing the best algorithm for given modulo and degree.
    Choice is based on heuristic measurments.

.. function:: void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f)

    Factorises a non-constant polynomial ``f`` into monic irreducible
    factors using the Cantor-Zassenhaus algorithm.

.. function:: void fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t poly)

    Factorises a non-constant polynomial ``poly`` into monic irreducible
    factors using the fast version of Cantor-Zassenhaus algorithm proposed by
    Kaltofen and Shoup (1998). More precisely this algorithm uses a
    baby step/giant step strategy for the distinct-degree factorization
    step. If ``flint_get_num_threads()`` is greater than one
    ``fmpz_mod_poly_factor_distinct_deg_threaded`` is used.

.. function:: void fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors, const fmpz_mod_poly_t f)

    Factorises a non-constant polynomial ``f`` into monic irreducible
    factors using the Berlekamp algorithm.

.. function:: void * _fmpz_mod_poly_interval_poly_worker(void* arg_ptr)

    Worker function to compute interval polynomials in distinct degree
    factorisation. Input/output is stored in
    ``fmpz_mod_poly_interval_poly_arg_t``.


