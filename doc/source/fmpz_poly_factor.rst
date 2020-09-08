.. _fmpz-poly-factor:

**fmpz_poly_factor.h** -- factorisation of polynomials over the integers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_poly_factor_struct

.. type:: fmpz_poly_factor_t

    Description.

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

    Assumes that `\deg(p) \geq 2` and `e \neq 0`.

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


.. function:: void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, fmpz_poly_t F)

    Takes as input a polynomial `F` and a freshly initialized factor 
    structure ``fac``.  Updates ``fac`` to contain a factorization 
    of `F` into (not necessarily irreducible) factors that themselves 
    have no repeated factors.  None of the returned factors will have 
    the same exponent. That is we return `g_i` and unique `e_i` such that 

    .. math ::


        F = c \prod_{i} g_i^{e_i}


    where `c` is the signed content of `F` and `\gcd(g_i, g_i') = 1`.

.. function:: void fmpz_poly_factor_zassenhaus_recombination(fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac, const fmpz_poly_t F, const fmpz_t P, slong exp)

    Takes as input a factor structure ``lifted_fac`` containing a 
    squarefree factorization of the polynomial `F \bmod p`. The algorithm 
    does a brute force search for irreducible factors of `F` over the 
    integers, and each factor is raised to the power ``exp``.

    The impact of the algorithm is to augment a factorization of 
    ``F^exp`` to the factor structure ``final_fac``.

.. function:: void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, slong exp, fmpz_poly_t f, slong cutoff, int use_van_hoeij)

    This is the internal wrapper of Zassenhaus.

    It will attempt to find a small prime such that `f` modulo `p` has 
    a minimal number of factors.  If it cannot find a prime giving less 
    than ``cutoff`` factors it aborts.  Then it decides a `p`-adic 
    precision to lift the factors to, hensel lifts, and finally calls 
    Zassenhaus recombination.

    Assumes that `\operatorname{len}(f) \geq 2`.

    Assumes that `f` is primitive.

    Assumes that the constant coefficient of `f` is non-zero.  Note that 
    this can be easily achieved by taking out factors of the form `x^k` 
    before calling this routine.

    If the final flag is set, the function will use the van Hoeij factorisation
    algorithm with gradual feeding and mod `2^k` data truncation to find
    factors when the number of local factors is large.

.. function:: void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, fmpz_poly_t F)

    A wrapper of the Zassenhaus factoring algorithm, which takes as input 
    any polynomial `F`, and stores a factorization in ``final_fac``.

    The complexity will be exponential in the number of local factors 
    we find for the components of a squarefree factorization of `F`.

.. function:: void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)
              void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)

    Inserts the factorisation of the quadratic (resp. cubic) polynomial *f* into *fac* with
    multiplicity *exp*. This function requires that the content of *f* has
    been removed, and does not update the content of *fac*.
    The factorzation is calculated over `\mathbb{R}` or `\mathbb{Q}_2` and then tested over `\mathbb{Z}`.

.. function:: void fmpz_poly_factor(fmpz_poly_factor_t final_fac, fmpz_poly_t F)

    A wrapper of the Zassenhaus and van Hoeij factoring algorithms, which takes
    as input any polynomial `F`, and stores a factorization in
    ``final_fac``.

