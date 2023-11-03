.. _fmpzi:

**fmpzi.h** -- Gaussian integers
===============================================================================

This module allows working with elements of the ring `\mathbb{Z}[i]`.
At present, only a minimal interface is provided.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpzi_struct

.. type:: fmpzi_t

    Contains a pairs of integers representing the real and imaginary
    parts. An *fmpzi_t* is defined as an array
    of length one of type *fmpzi_struct*, permitting an *fmpzi_t* to
    be passed by reference.

.. macro:: fmpzi_realref(x)

    Macro giving a pointer to the real part of *x*.

.. macro:: fmpzi_imagref(x)

    Macro giving a pointer to the imaginary part of *x*.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void fmpzi_init(fmpzi_t x)

.. function:: void fmpzi_clear(fmpzi_t x)

.. function:: void fmpzi_swap(fmpzi_t x, fmpzi_t y)

.. function:: void fmpzi_zero(fmpzi_t x)

.. function:: void fmpzi_one(fmpzi_t x)

.. function:: void fmpzi_set(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_set_si_si(fmpzi_t res, slong a, slong b)

Input and output
-------------------------------------------------------------------------------

.. function:: void fmpzi_print(const fmpzi_t x)

Random number generation
-------------------------------------------------------------------------------

.. function:: void fmpzi_randtest(fmpzi_t res, flint_rand_t state, mp_bitcnt_t bits)

Properties
-------------------------------------------------------------------------------

.. function:: int fmpzi_equal(const fmpzi_t x, const fmpzi_t y)

.. function:: int fmpzi_is_zero(const fmpzi_t x)

.. function:: int fmpzi_is_one(const fmpzi_t x)

Units
-------------------------------------------------------------------------------

.. function:: int fmpzi_is_unit(const fmpzi_t x)

.. function:: slong fmpzi_canonical_unit_i_pow(const fmpzi_t x)

.. function:: void fmpzi_canonicalise_unit(fmpzi_t res, const fmpzi_t x)

Norms
-------------------------------------------------------------------------------

.. function:: slong fmpzi_bits(const fmpzi_t x)

.. function:: void fmpzi_norm(fmpz_t res, const fmpzi_t x)

Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpzi_conj(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_neg(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sqr(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)

Division
-------------------------------------------------------------------------------

.. function:: void fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)

    Sets *q* to the quotient of *x* and *y*, assuming that the
    division is exact.

.. function:: void fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

    Computes a quotient and remainder satisfying
    `x = q y + r` with `N(r) \le N(y)/2`, with a canonical
    choice of remainder when breaking ties.

.. function:: void fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

    Computes a quotient and remainder satisfying
    `x = q y + r` with `N(r) < N(y)`, with an implementation-defined,
    non-canonical choice of remainder.

.. function:: slong fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x)

    Divide *x* exactly by the largest possible power `(1+i)^k`
    and return the exponent *k*.

GCD
-------------------------------------------------------------------------------

.. function:: void fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_shortest(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

    Computes the GCD of *x* and *y*. The result is in canonical
    unit form.

    The *euclidean* version is a straightforward implementation
    of Euclid's algorithm. The *euclidean_improved* version is
    optimized by performing approximate divisions.
    The *binary* version uses a (1+i)-ary analog of the binary
    GCD algorithm for integers [Wei2000]_.
    The *shortest* version finds the GCD as the shortest vector in a lattice.
    The default version chooses an algorithm automatically.

Primality testing
--------------------------------------------------------------------------------

.. function:: int fmpzi_is_prime(const fmpzi_t n)

    Check whether `n` is a Gaussian prime.

.. function:: int fmpzi_is_probabprime(const fmpzi_t n)

    Check whether `n` is a probable Gaussian prime.
