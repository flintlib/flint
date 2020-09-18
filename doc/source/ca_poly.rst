.. _ca-poly:

**ca_poly.h** -- dense univariate polynomials over the real and complex numbers
===============================================================================

An :type:`ca_poly_t` represents a dense univariate
polynomial over the real or complex numbers,
implemented as an array of coefficients of type :type:`ca_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Warnings:

* A polynomial is always normalised by removing leading zero coefficients.
  Coefficients will not be removed when Calcium is unable to prove
  that they are zero. The represented degree can therefore be larger
  than the degree of the mathematical polynomial.
  When the correct degree is needed, it is important to verify
  the leading coefficient.

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

    Strips any trailing coefficients which are identical to zero.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx)

    Prints *poly* to standard output. The coefficients are printed on separate lines.

.. function:: void ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx)

    Prints a decimal representation of *poly* with precision specified by *digits*.
    The coefficients are comma-separated and the whole list is enclosed in square brackets.

.. raw:: latex

    \newpage
