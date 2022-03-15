.. generic-rings documentation master file, created by
   sphinx-quickstart on Sun Mar 13 11:37:54 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

generic-rings documentation
=========================================

This module provides a plain C implementation of generic rings
intended for use with Flint, Antic, Arb and Calcium.
The idea is to support fully generic recursive constructions of rings
(polynomials, fractions, power series, matrices, etc.) over arbitrary
user-defined base rings (including all Flint/Antic/Arb/Calcium types)
with the ability to add specializations (e.g. ``fmpz_mat`` functions
instead of generic functions for matrices over ``fmpz``).
This code is *experimental*.

.. contents:: Table of contents
   :depth: 2
   :local:
   :backlinks: none


Design goals
-------------------------------------------------------------------------------

We use element pointers + context objects holding method pointers.
Vectors of elements can be packed densely without indirection.

Principles/goals/benefits
...............................................................................

* Plain C, similar interface to existing Flint code.
* Small code size, fast compilation.
* Possible to pack data efficiently (down to 1 byte / element).
* Data layouts backwards compatible with most existing Flint element types,
  and mostly also with polynomials, matrices, etc.
* Support all unusual cases in Flint/Arb/Calcium uniformly (error handling,
  inexact rings, noncomputable rings, context objects), with a uniform
  interface.
* Support fast stack based allocation of temporary variables and arrays
  (with safe limits).
* Fully in-place operations.
* Very fast (a few cycles) runtime construction of rings / context objects.
  (May be hundreds of cycles when constructing a new
  method table, but this is rarely needed on the fly.)
* Possibility to use generic default methods or provide optimized versions
  (e.g. for vector operations).

Disadvantages
...............................................................................

* Runtime generics in C means very little compiler protection against
  type errors; rigorous testing is crucial.
* Some function signatures become more complicated in order to provide
  uniform error handling.
* A few cycles overhead for each dynamic function lookup and
  function call. Function pointer generics also preclude intra-method
  compiler optimizations and inlining.
  Overloading vector methods should partially compensate.

Possible applications
...............................................................................

* At minimum, this could be useful to bootstrap new types: we only
  need to provide some basic methods (``init``, ``clear``, ``set_fmpz``,
  ``randtest``, ``add``, ``mul``, ``equal``, etc.) and the generic
  code provides fallbacks for everything else (which can be overloaded
  later for better performance). Importantly, the generic code also
  provides unit tests including tests for boring and easily borked
  things like type conversions.
* This could simplify interfacing from other libraries and languages;
  it is enough to wrap the generic interface once, without having
  to wrap all the methods for every single Flint type (with all
  their quirks and subtle differences in interface).
* There is already plenty of ad-hoc generic code in Flint
  and its descendants:

  * Function pointer generics are used, for example, in some ``mpoly`` code and in ``fmpz_mod``.

  * Template-based generics for the different ``fq`` implementations.

  * Type union + switch statement generics for ``fq_default`` and ``nf``.

  * Some Arb functions use ugly hacks to distinguish between constants and power series, etc.

  The generics module could potentially replace such code and also
  eliminate a lot of other copy-pasted boilerplate.
  By specializing methods at runtime for rings with different parameters,
  it should even be possible to improve performance in some cases.


Error handling
-------------------------------------------------------------------------------

To implement computations over a ring `R`,
it is sometimes convenient to extend the ring to a set
`R' = R \cup \{ \text{undefined}, \text{unknown} \}`.
An *undefined* (error) value allows us to extend partial functions
to total functions.
Alternatively,
we could use some arbitrary default value in `R`,
say `\text{undefined} = 0`; this is often done in
formal theorem provers,
but it may be undesirable in a regular programming
environment as it makes it harder to detect bugs.
An *unknown* value is useful in cases where a result
may exist in principle but cannot be computed.

Ideally, we would represent `R'` as a type-level extension of `R`,
but this tricky in C since we would either have to
wrap elements in a larger structure
or reserve bit patterns in each type for special values.
In any case, it is useful to assume in low-level code
that ring elements really represent ring elements
so that there are fewer special cases to handle.
We also need some form of error handling for conversions
to standard C types.
For these reasons, we handle special values (undefined, unknown)
using return codes.

Functions can return a combination of the following status flags:

.. macro:: GR_SUCCESS

    The operation finished as expected, i.e. the result
    is a correct element of the target type.

.. macro:: GR_DOMAIN

    The result does not have a value in the domain of the target
    ring or type, i.e. the result is mathematically undefined.
    This occurs, for example, on division by zero
    or when attempting to compute the square root of a non-square.
    It also occurs when attempting to convert a too large value
    to a bounded type (example: ``get_ui()``
    with input `n \ge 2^{64}`).

.. macro:: GR_UNABLE

    The operation could not be performed because
    of limitations of the implementation or the data representation,
    i.e. the result is unknown. Typical reasons:

    * The result would be too large to fit in memory
    * The inputs are inexact and an exact comparison is needed
    * The computation would take too long
    * An algorithm is not yet implemented for this case

    If this flag is set, there is also potentially a domain error
    (but this is unknown).

.. macro:: GR_WRONG

    Test failure. This is only used in test code.

When the status code is any other value than ``GR_SUCCESS``, any
output variables may be set to meaningless values.

For uniformity, even functions that should never fail return a status
code (we might want to wrap such functions in asserts).
Flags can be OR'ed and checked only at the top level of a computation
to avoid complex control flow.

Main types
-------------------------------------------------------------------------------

.. type:: gr_ptr

    Pointer to a ring element. This is an alias for ``void *``
    so that it can be used with any C type.

.. type:: gr_srcptr

    Pointer to a read-only ring element. This is an alias for
    ``const void *`` so that it can be used with any C type.

.. type:: gr_ctx_struct

.. type:: gr_ctx_t

    A context object representing a mathematical ring *R*.
    It contains the following data:

    * Flags describing useful properties of the ring.
    * The size (number of bytes) of each element.
    * A pointer to a method table.
    * Optionally a pointer to data defining parameters of the ring
      (e.g. modulus of a residue ring; element ring and dimensions
      of a matrix ring; precision of an inexact ring).

    A :type:`gr_ctx_t` is defined as an array of length one of type
    :type:`gr_ctx_struct`, permitting a :type:`gr_ctx_t` to be
    passed by reference.
    Context objects are not normally passed as ``const`` in order
    to allow storing mutable caches, additional
    debugging information, etc.

.. type:: gr_ctx_ptr

    Pointer to a context object.

Observe that there is no type to represent a single generic element
as a struct since we do not know the size of a generic element at
compile time.
Memory for single elements can either be allocated on the stack
with the special macros provided below, or as usual with ``malloc``.

When using generic methods with a known type like
``fmpz_t``, the usual type can of course be used.
Users may wish to define their own union types when only some
particular types will appear in an application.

Ring constructions
-------------------------------------------------------------------------------

Builtin base rings
...............................................................................

.. function:: void gr_ctx_init_fmpq(gr_ctx_t ctx)

    Initializes *ctx* to the field of rational numbers
    `\mathbb{Q}` with elements of type :type:`fmpq`.

.. function:: void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`uint8`. We require `1 \le n \le 255`.

.. function:: void gr_ctx_init_real_qqbar(gr_ctx_t ctx)
              void gr_ctx_init_complex_qqbar(gr_ctx_t ctx)

    Initializes *ctx* to the field of real or complex algebraic
    numbers with elements of type :type:`qqbar`.


Derived rings
...............................................................................

.. function:: void gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

    Initializes *ctx* to the ring of densely stored *n* by *n* matrices
    over the given *base_ring*.
    Elements have type :type:`gr_mat_struct`.

Operations
-------------------------------------------------------------------------------

.. function:: int gr_ctx_clear(gr_ctx_t ctx)

    Clears the ring context object *ctx*, freeing any memory
    allocated by this object.

    Some rings may require that no elements are cleared after calling
    this method, and may leak memory if not all elements have
    been cleared when calling this method.

    If *ctx* is derived from a base ring, the base ring context
    may also be required to stay alive until after this
    method is called.

.. function:: int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx)

    Writes a description of the ring *ctx* to the stream *out*.

.. function:: int gr_init(gr_ptr res, gr_ctx_t ctx)

    Initializes *res* to a valid variable and sets it to the
    zero element of the ring *ctx*.

.. function:: int gr_clear(gr_ptr res, gr_ctx_t ctx)

    Clears *res*, freeing any memory allocated by this object.

.. function:: int gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)

    Swaps *x* and *y* efficiently.

.. function:: int gr_randtest(gr_ptr res, flint_rand_t state, const void * options, gr_ctx_t ctx)

    Sets *res* to a random element of the ring.

.. function:: int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)

    Writes a representation of *x* to the stream *out*.

.. function:: int gr_zero(gr_ptr res, gr_ctx_t ctx)
              int gr_one(gr_ptr res, gr_ctx_t ctx)
              int gr_neg_one(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to the element 0, 1 or -1 of the ring.

.. function:: int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a copy of the element *x*.

.. function:: int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx)
              int gr_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
              int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
              int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx)

    Sets *res* to the image of the integer or rational number *x*
    in the ring *ctx*.
    The *fmpq* method may return the flag ``GR_DOMAIN`` if the
    denominator of *x* is not invertible.

.. function:: int gr_is_zero(int * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_is_one(int * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_is_neg_one(int * res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to 1 if *x* is equal to the element 0, 1 or -1 of the
    ring, respectively, and sets *res* to 0 otherwise.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the comparison.

.. function:: int gr_equal(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to 1 if the elements *x* and *y* are equal,
    and sets *res* to 0 otherwise.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the comparison.

Arithmetic
........................................................................

.. function:: int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

.. function:: int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

.. function:: int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

Division
........................................................................

The default implementations of the following methods check for divisors
0, 1, -1 and otherwise return ``GR_UNABLE``.
Particular rings should override the methods when an inversion
or division algorithm is available.
The base rings corresponding to
the following types have complete algorithms
to detect inverses and compute quotients: ``fmpq``, ``qqbar``, ``nmod8``.

.. function:: int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

.. function:: int gr_is_invertible(int * res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to 1 if *x* has a multiplicative inverse in the present ring
    (i.e. if *x* is a unit),
    and sets *res* to 0 if *x* does not have a multiplicative inverse.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the computation.

.. function:: int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the multiplicative inverse of *x* in the present ring,
    if such an element exists.
    Returns the flag ``GR_DOMAIN`` if *x* is not invertible, or
    ``GR_UNABLE`` if the implementation is unable to perform
    the computation.

Powering
........................................................................

.. function:: int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to the power `x ^ y`.
    Returns the flag ``GR_DOMAIN`` if this power cannot be assigned
    a meaningful value in the present ring, or ``GR_UNABLE`` if
    the implementation is unable to perform the computation.

    Default implementations of the powering methods support raising
    elements to integer powers using a generic implementation of
    exponentiation by squaring. Particular rings
    should override these methods with faster versions or
    to support more general notions of exponentiation when possible.

Square roots
........................................................................

The default implementations of the following methods check for the
elements 0 and 1 and otherwise return ``GR_UNABLE``.
Particular rings should override the methods when a square
root algorithm is available.
The base rings corresponding to
the following types have complete algorithms
to detect squares and compute square roots: ``fmpq``, ``qqbar``.

In subrings of `\mathbb{C}`, it is implied that the principal
square root is computed; in other cases (e.g. in finite fields),
the choice of root is implementation-dependent.

.. function:: int gr_is_square(int * res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to 1 if *x* is a perfect square in the present ring,
    and sets *res* to 0 if *x* it not a perfect square.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the computation.

.. function:: int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a square root of *x* (respectively reciprocal
    square root) in the present ring, if such an element exists.
    Returns the flag ``GR_DOMAIN`` if *x* is not a perfect square
    (also for zero, when computing the reciprocal square root), or
    ``GR_UNABLE`` if the implementation is unable to perform
    the computation.
