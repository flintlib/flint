.. _fmpz-poly-factor:

**fmpz_poly_factor.h** -- factorisation of polynomials over the integers
===============================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_poly_factor_struct

.. type:: fmpz_poly_factor_t

Memory management
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_factor_init(fmpz_poly_factor_t fac)

    Initialises a new factor structure.

.. function:: void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, slong alloc)

    Initialises a new factor structure, providing space for
    at least ``alloc`` factors.

.. function:: void fmpz_poly_factor_realloc(fmpz_poly_factor_t fac, slong alloc)

    Reallocates the factor structure to provide space for
    precisely ``alloc`` factors.

.. function:: void fmpz_poly_factor_fit_length(fmpz_poly_factor_t fac, slong len)

    Ensures that the factor structure has space for at
    least ``len`` factors.  This functions takes care
    of the case of repeated calls by always at least
    doubling the number of factors the structure can hold.

.. function:: void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)

    Releases all memory occupied by the factor structure.


Manipulating factors
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_factor_set(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)

    Sets ``res`` to the same factorisation as ``fac``.

.. function:: void fmpz_poly_factor_insert(fmpz_poly_factor_t fac, const fmpz_poly_t p, slong e)

    Adds the primitive polynomial `p^e` to the factorisation ``fac``.

    Assumes that `\deg(p) \geq 2` and `e \neq 0`. These conditions are not
    checked.

.. function:: void fmpz_poly_factor_concat(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)

    Concatenates two factorisations.

    This is equivalent to calling :func:`fmpz_poly_factor_insert`
    repeatedly with the individual factors of ``fac``.

    Does not support aliasing between ``res`` and ``fac``.


Input and output
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_factor_print(const fmpz_poly_factor_t fac)

    Prints the entries of ``fac`` to standard output.


Factoring algorithms
--------------------------------------------------------------------------------


.. function:: void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, const fmpz_poly_t F)

    Takes as input a polynomial `F` and a freshly initialized factor
    structure ``fac``.  Updates ``fac`` to contain a factorization
    of `F` into (not necessarily irreducible) factors that themselves
    have no repeated factors.  None of the returned factors will have
    the same exponent. That is we return `g_i` and unique `e_i` such that

    .. math::


        F = c \prod_{i} g_i^{e_i}


    where `c` is the signed content of `F` and `\gcd(g_i, g_i') = 1`.

.. function:: void fmpz_poly_factor_zassenhaus_recombination(fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac, const fmpz_poly_t F, const fmpz_t P, slong exp)

    Takes as input a factor structure ``lifted_fac`` containing a
    squarefree factorization of the polynomial `F \bmod p`. The algorithm
    does a brute force search for irreducible factors of `F` over the
    integers, and each factor is raised to the power ``exp``.

    The impact of the algorithm is to augment a factorization of
    ``F^exp`` to the factor structure ``final_fac``.

.. function:: void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij)

    This is the internal wrapper of Zassenhaus.

    It will attempt to find a small prime such that `f` modulo `p` has
    a minimal number of factors.  If it cannot find a prime giving less
    than ``cutoff`` factors it aborts.  Then it decides a `p`-adic
    precision to lift the factors to, Hensel lifts, and finally calls
    Zassenhaus recombination.

    Assumes that `\operatorname{len}(f) \geq 2`.

    Assumes that `f` is primitive.

    Assumes that the constant coefficient of `f` is non-zero.  Note that
    this can be easily achieved by taking out factors of the form `x^k`
    before calling this routine.  Neither the primitivity of `f` nor the
    non-zero constant coefficient is checked.

    If the final flag is set, the function will use the van Hoeij factorisation
    algorithm with gradual feeding and mod `2^k` data truncation to find
    factors when the number of local factors is large.

.. function:: void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, const fmpz_poly_t F)

    A wrapper of the Zassenhaus factoring algorithm, which takes as input
    any polynomial `F`, and stores a factorization in ``final_fac``.

    The complexity will be exponential in the number of local factors
    we find for the components of a squarefree factorization of `F`.

.. function:: void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)
              void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)

    Inserts the factorisation of the quadratic (resp. cubic) polynomial *f* into *fac* with
    multiplicity *exp*. This function requires that the content of *f* has
    been removed (this is not checked), and does not update the content of
    *fac*.
    The factorization is calculated over `\mathbb{R}` or `\mathbb{Q}_2` and then tested over `\mathbb{Z}`.

.. function:: void fmpz_poly_factor(fmpz_poly_factor_t final_fac, const fmpz_poly_t F)

    A wrapper of the Zassenhaus and van Hoeij factoring algorithms, which takes
    as input any polynomial `F`, and stores a factorization in
    ``final_fac``.

.. function:: int _fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong p)

    Given `T` irreducible over `\mathbb{Q}` of degree at least 1 and a prime
    `p`, attempts to certify that the inflation `T(x^p)` is irreducible over
    `\mathbb{Q}`.  Returns 1 if a certificate was found, in which case
    `T(x^p)` is provably irreducible, and 0 if the search was inconclusive,
    in which case nothing is claimed (in particular, `T(x^p)` may still be
    irreducible).  The behaviour is undefined if `T` is reducible.  `T` need
    not be primitive and may have negative leading coefficient.

    By Capelli's theorem, for prime `p` the polynomial `T(x^p)` is
    irreducible over `\mathbb{Q}` if and only if `\theta`, the image of `x`
    in `K = \mathbb{Q}[x]/(T)`, is not a `p`-th power in `K`.  For
    `\deg T = 1`, and for `\deg T = 2` with `p = 2`, this condition is
    decided exactly (so the return value is 1 precisely when `T(x^p)` is
    irreducible).  Otherwise the function searches for a rational prime `q`,
    with `T` squarefree modulo `q` and `q` dividing neither the leading nor
    the constant coefficient of `T`, such that `\theta` reduces to a
    non-`p`-th power in some residue field of `K` above `q`; any such
    witness proves `\theta \notin K^p`.  The search is heuristically tuned:
    its budget is linear in `\deg T` (for inputs with abelian splitting
    fields of exponent `p`, such as the deflations of Swinnerton-Dyer
    polynomials, witnesses have density only `1/(p \deg T)` among rational
    primes and occur only in degree-1 residue fields), and it terminates
    early when accumulated local evidence makes irreducibility unlikely.

.. function:: int fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong d)

    Given `T` irreducible over `\mathbb{Q}` of degree at least 1 and
    `d \geq 1`, returns 1 if the inflation `T(x^d)` is certified irreducible
    over `\mathbb{Q}`, and 0 if the certification was inconclusive, in which
    case nothing is claimed.  Chains the prime-step certificate over the
    prime factorisation of `d`, via `T(x^d) = U(x^{d/p})` with
    `U = T(x^p)`.  The behaviour is undefined if `T` is reducible.
    This certificate is used by :func:`fmpz_poly_factor` to avoid
    refactoring inflations of irreducible deflation factors.

