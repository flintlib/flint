.. _ca-poly:

**ca_poly.h** -- dense univariate polynomials over the real and complex numbers
===============================================================================

A :type:`ca_poly_t` represents a univariate
polynomial over the real or complex numbers,
implemented as an array of coefficients of type :type:`ca_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Warnings:

* A polynomial is always normalised by removing zero coefficients
  at the top.
  Coefficients will not be removed when Calcium is unable to prove
  that they are zero. The represented degree can therefore be larger
  than the degree of the mathematical polynomial.
  When the correct degree is needed, it is important to verify
  the leading coefficient.
  (Of course, this will never be an issue with polynomials
  that are explicitly monic, for example.)

* The special values *Undefined*, unsigned infinity and signed infinity
  supported by the scalar :type:`ca_t` type
  are not really meaningful as coefficients of polynomials.
  We normally assume that the user does not assign those values to
  coefficients of polynomials, and the functions in this module will
  likewise normally not generate such coefficients.
  *Unknown* can still appear as a coefficient representing a
  number that is inaccessible for computation.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ca_poly_struct

.. type:: ca_poly_t

    Contains a pointer to an array of coefficients (*coeffs*), the used
    length (*length*), and the allocated size of the array (*alloc*).

    A *ca_poly_t* is defined as an array of length one of type
    *ca_poly_struct*, permitting an *ca_poly_t* to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void ca_poly_init(ca_poly_t poly, ca_ctx_t ctx)

    Initializes the polynomial for use, setting it to the zero polynomial.

.. function:: void ca_poly_clear(ca_poly_t poly, ca_ctx_t ctx)

    Clears the polynomial, deallocating all coefficients and the
    coefficient array.

.. function:: void ca_poly_fit_length(ca_poly_t poly, slong len, ca_ctx_t ctx)

    Makes sure that the coefficient array of the polynomial contains at
    least *len* initialized coefficients.

.. function:: void _ca_poly_set_length(ca_poly_t poly, slong len, ca_ctx_t ctx)

    Directly changes the length of the polynomial, without allocating or
    deallocating coefficients. The value should not exceed the allocation length.

.. function:: void _ca_poly_normalise(ca_poly_t poly, ca_ctx_t ctx)

    Strips any top coefficients which can be proved identical to zero.

Assignment and simple values
-------------------------------------------------------------------------------

.. function:: void ca_poly_zero(ca_poly_t poly, ca_ctx_t ctx)

    Sets *poly* to the zero polynomial.

.. function:: void ca_poly_one(ca_poly_t poly, ca_ctx_t ctx)

    Sets *poly* to the constant polynomial 1.

.. function:: void ca_poly_x(ca_poly_t poly, ca_ctx_t ctx)

    Sets *poly* to the monomial *x*.

.. function:: void ca_poly_set_ca(ca_poly_t poly, const ca_t c, ca_ctx_t ctx)

.. function:: void ca_poly_set_si(ca_poly_t poly, slong c, ca_ctx_t ctx)

    Sets *poly* to the constant polynomial *c*.

.. function:: void ca_poly_set(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)

.. function:: void ca_poly_set_fmpz_poly(ca_poly_t res, const fmpz_poly_t src, ca_ctx_t ctx)

.. function:: void ca_poly_set_fmpq_poly(ca_poly_t res, const fmpq_poly_t src, ca_ctx_t ctx)

    Sets *poly* the polynomial *src*.

Random generation
-------------------------------------------------------------------------------

.. function:: void ca_poly_randtest(ca_poly_t poly, flint_rand_t state, slong len, slong depth, slong bits, ca_ctx_t ctx)

    Sets *poly* to a random polynomial of length up to *len* and with entries having complexity up to
    *depth* and *bits* (see :func:`ca_randtest`).

.. function:: void ca_poly_randtest_rational(ca_poly_t poly, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)

    Sets *poly* to a random rational polynomial of length up to *len* and with entries up to *bits* bits in size.


Input and output
-------------------------------------------------------------------------------

.. function:: void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx)

    Prints *poly* to standard output. The coefficients are printed on separate lines.

.. function:: void ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx)

    Prints a decimal representation of *poly* with precision specified by *digits*.
    The coefficients are comma-separated and the whole list is enclosed in square brackets.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_poly_neg(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)

    Sets *res* to the negation of *src*.

.. function:: void _ca_poly_add(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)

.. function:: void ca_poly_add(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

    Sets *res* to the sum of *poly1* and *poly2*.

.. function:: void _ca_poly_sub(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)

.. function:: void ca_poly_sub(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

    Sets *res* to the difference of *poly1* and *poly2*.

.. function:: void _ca_poly_mul(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)

.. function:: void ca_poly_mul(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

    Sets *res* to the product of *poly1* and *poly2*.

.. function:: void ca_poly_mul_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)

    Sets *res* to *poly1* multiplied by the scalar *c*.

.. function:: void _ca_poly_mullow(ca_ptr C, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, slong n, ca_ctx_t ctx)

.. function:: void ca_poly_mullow(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, slong n, ca_ctx_t ctx)

    Sets *res* to the product of *poly1* and *poly2* truncated to length *n*.

Roots and factorization
-------------------------------------------------------------------------------

.. function:: void _ca_poly_set_roots(ca_ptr poly, ca_srcptr roots, slong n, ca_ctx_t ctx)

.. function:: void ca_poly_set_roots(ca_poly_t poly, ca_vec_t roots, ca_ctx_t ctx)

    Sets *poly* to the monic polynomial with the *n* roots
    given in the vector *roots*. That is, sets *poly* to
    `(x-r_0)(x-r_1)\cdots(x-r_{n-1})`.

.. function:: int _ca_poly_roots(ca_ptr roots, ca_srcptr poly, slong len, ca_ctx_t ctx)

.. function:: int ca_poly_roots(ca_vec_t roots, const ca_poly_t poly, ca_ctx_t ctx)

    Attempts to compute all complex roots of the given polynomial *poly*.
    On success, returns 1 and sets *roots* to the vector of roots.
    On failure, returns 0 and leaves the values in *roots* arbitrary.

    The roots are returned in arbitrary order, but repeated according
    to their multiplicity.

    Failure will occur if the leading coefficient of *poly* cannot
    be proved to be nonzero, if determining the correct multiplicities
    fails, or if the builtin algorithms do not have a means to
    represent the roots symbolically.



.. raw:: latex

    \newpage
