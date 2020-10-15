.. _ca-vec:

**ca_vec.h** -- vectors of real and complex numbers
===============================================================================

A :type:`ca_vec_t` represents a vector of real or complex numbers,
implemented as an array of coefficients of type :type:`ca_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients
(taking :type:`ca_ptr` and :type:`ca_srcptr` arguments),
and a non-underscore method which takes :type:`ca_vec_t` input
and performs automatic memory management.

Unlike :type:`ca_poly_t`, a :type:`ca_vec_t` is not normalised
by removing zero coefficients; it retains the exact length
assigned by the user.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ca_vec_struct

.. type:: ca_vec_t

    Contains a pointer to an array of entries (*coeffs*), the used
    length (*length*), and the allocated size of the array (*alloc*).

    A *ca_vec_t* is defined as an array of length one of type
    *ca_vec_struct*, permitting an *ca_vec_t* to
    be passed by reference.

.. macro:: ca_vec_entry(vec, i)

    Macro returning a pointer to entry *i* in the vector *vec*.
    The index must be in bounds.

Memory management
-------------------------------------------------------------------------------

.. function:: ca_ptr _ca_vec_init(slong len, ca_ctx_t ctx)

    Returns a pointer to an array of *len* coefficients
    initialized to zero.

.. function:: void ca_vec_init(ca_vec_t vec, slong len, ca_ctx_t ctx)

    Initializes *vec* to a length *len* vector. All entries
    are set to zero.

.. function:: void _ca_vec_clear(ca_ptr vec, slong len, ca_ctx_t ctx)

    Clears all *len* entries in *vec* and frees the pointer
    *vec* itself.

.. function:: void ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx)

    Clears the vector *vec*.

.. function:: void ca_vec_swap(ca_vec_t vec1, ca_vec_t vec2, ca_ctx_t ctx)

    Swaps the vectors *vec1* and *vec2* efficiently.

Length
-------------------------------------------------------------------------------

.. function:: slong ca_vec_length(const ca_vec_t vec, ca_ctx_t ctx)

    Returns the length of *vec*.

.. function:: void ca_vec_set_length(ca_vec_t vec, slong len, ca_ctx_t ctx)

    Sets the length of *vec* to *len*.
    If *vec* is shorter on input, it will be zero-extended.
    If *vec* is longer on input, it will be truncated.

Assignment
-------------------------------------------------------------------------------

.. function:: void _ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)

    Sets *res* to a copy of *src* of length *len*.

.. function:: void ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)

    Sets *res* to a copy of *src*.

Special vectors
-------------------------------------------------------------------------------

.. function:: void _ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx)

    Sets the *len* entries in *res* to zeros.

.. function:: void ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx)

    Sets *res* to the length *len* zero vector.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_vec_print(const ca_vec_t vec, ca_ctx_t ctx)

    Prints *vec* to standard output. The coefficients are printed on separate lines.

.. function:: void ca_vec_printn(const ca_vec_t poly, slong digits, ca_ctx_t ctx)

    Prints a decimal representation of *vec* with precision specified by *digits*.
    The coefficients are comma-separated and the whole list is enclosed in square brackets.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void _ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)

.. function:: void ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)

    Sets *res* to the negation of *src*.

.. function:: void _ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)

.. function:: void _ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)

    Sets *res* to the sum or difference of *vec1* and *vec2*,
    all vectors having length *len*.

.. function:: void _ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

    Sets *res* to *src* multiplied by *c*, all vectors having
    length *len*.

.. function:: void _ca_vec_scalar_div_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

    Sets *res* to *src* divided by *c*, all vectors having
    length *len*.

.. function:: void _ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

    Adds *src* multiplied by *c* to the vector *res*, all vectors having
    length *len*.

.. function:: void _ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

    Subtracts *src* multiplied by *c* from the vector *res*, all vectors having
    length *len*.

Internal representation
---------------------------------------------------------------------------------

.. function:: int _ca_vec_is_fmpq_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)

    Checks if all elements of *vec* are structurally rational numbers.

.. function:: int _ca_vec_fmpq_vec_is_fmpz_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)

    Assuming that all elements of *vec* are structurally rational numbers,
    checks if all elements are integers.

.. function:: void _ca_vec_fmpq_vec_get_fmpz_vec_den(fmpz * c, fmpz_t den, ca_srcptr vec, slong len, ca_ctx_t ctx)

    Assuming that all elements of *vec* are structurally rational numbers,
    converts them to a vector of integers *c* on a common denominator
    *den*.

.. function:: void _ca_vec_set_fmpz_vec_div_fmpz(ca_ptr res, const fmpz * v, const fmpz_t den, slong len, ca_ctx_t ctx)

    Sets *res* to the rational vector given by numerators *v*
    and the common denominator *den*.

.. raw:: latex

    \newpage
