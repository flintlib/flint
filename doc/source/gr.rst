.. _gr:

**gr.h** -- generic rings
===============================================================================

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

For uniformity, most operations return a status code, even operations
that are not typically expected to fail (we might want to wrap
such functions in asserts).

* Pure "container" operations like ``init``, ``clear`` and ``swap``
  do not return a status code.

* Pure predicate functions (see below)
  return ``T_TRUE`` / ``T_FALSE`` / ``T_UNKNOWN``
  instead of computing a separate boolean value and error code.

Flags can be OR'ed and checked only at the top level of a computation
to avoid complex control flow.

.. macro:: MUST_SUCCEED(expr)

    Evaluates *expr* and asserts that the return value is
    ``GR_SUCCESS``.

Ring predicates
-------------------------------------------------------------------------------

We use the following type (borrowed from Calcium) instead of a C int
to represent boolean results, allowing the possibility
that the value is not computable:

.. enum:: truth_t

    Represents one of the following truth values:

    .. macro:: T_TRUE

    .. macro:: T_FALSE

    .. macro:: T_UNKNOWN

    Warning: the constants ``T_TRUE`` and ``T_FALSE`` do not correspond to 1 and 0.
    It is erroneous to write, for example ``!t`` if ``t`` is a 
    :type:`truth_t`. One should instead write ``t != T_TRUE``, ``t == T_FALSE``,
    etc. depending on whether the unknown case should be included
    or excluded.


Main types
-------------------------------------------------------------------------------

.. type:: gr_ptr

    Pointer to a ring element or array of contiguous ring elements.
    This is an alias for ``void *`` so that it can be used with any
    C type.

.. type:: gr_srcptr

    Pointer to a read-only ring element or read-only array of
    contiguous ring elements. This is an alias for
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

Base rings
...............................................................................

.. function:: void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)

    Initializes *ctx* to a random ring. This will currently
    only generate base rings.

.. function:: void gr_ctx_init_fmpz(gr_ctx_t ctx)

    Initializes *ctx* to the ring of integers
    `\mathbb{Z}` with elements of type :type:`fmpz`.

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
    numbers with elements of type :type:`qqbar_t`.

.. function:: void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
              void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the field of real or complex
    numbers represented by balls or boxes of type :type:`arb_t`
    and  :type:`acb_t`.

.. function:: void gr_ctx_init_real_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_ca(gr_ctx_t ctx)
              void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)

    Initializes *ctx* to the field of real, complex, real algebraic
    or complex algebraic numbers represented by elements of type
    :type:`ca_t`.


Derived rings
...............................................................................

.. function:: void gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

    Initializes *ctx* to a ring of densely stored *n* by *n* matrices
    over the given *base_ring*.
    Elements have type :type:`gr_mat_struct`.

.. function:: void gr_ctx_init_polynomial(gr_ctx_t ctx, gr_ctx_t base_ring)

    Initializes *ctx* to a ring of densely stored univariate polynomials
    over the given *base_ring*.
    Elements have type :type:`gr_poly_struct`.

Context operations
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
              int gr_ctx_print(gr_ctx_t ctx)
              int gr_ctx_println(gr_ctx_t ctx)
              int gr_ctx_get_str(char ** s, gr_ctx_t ctx)

    Writes a description of the ring *ctx* to the stream *out*,
    prints it to *stdout*, or sets *s* to a pointer to
    a heap-allocated string of the description (the user must free
    the string with ``flint_free``).
    The *println* version prints a trailing newline.

Element operations
--------------------------------------------------------------------------------

Memory management
................................................................................

.. function:: void gr_init(gr_ptr res, gr_ctx_t ctx)

    Initializes *res* to a valid variable and sets it to the
    zero element of the ring *ctx*.

.. function:: void gr_clear(gr_ptr res, gr_ctx_t ctx)

    Clears *res*, freeing any memory allocated by this object.

.. function:: void gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)

    Swaps *x* and *y* efficiently.

The following macros support allocating temporary variables efficiently.
Data will be allocated on the stack using ``alloca`` unless
the size is excessive (risking stack overflow), in which case
the implementation transparently switches to ``malloc``/``free``
instead. The usage pattern is as follows::

    {
        gr_ptr x, y;
        GR_TMP_START;

        GR_TMP_INIT2(x1, x2, ctx);

        /* do computations with x1, x2 */

        GR_TMP_CLEAR2(x1, x2, ctx);
        GR_TMP_END;
    }

Temporary allocations must be enclosed by the ``GR_TMP_START`` and
``GR_TMP_END`` markers, which should only occur
once in a block. In between, there
may be multiple calls to different init macros with matching clear
macros.
*Warning:* never use these macros directly inside a loop.
This is likely to overflow the stack, as memory will not
be reclaimed until the function exits.
Instead, allocate the needed space before entering
any loops, move the loop body to a separate function,
or allocate the memory on the heap if needed.

.. macro:: GR_TMP_START
           GR_TMP_END

    Markers for a block of temporary allocations.

.. macro:: GR_TMP_INIT_VEC(vec, len, ctx)
           GR_TMP_CLEAR_VEC(vec, len, ctx)

    Allocates and frees a vector of *len* contiguous elements, all
    initialized to the value 0, assigning the first element
    to the pointer *vec*.

.. macro:: GR_TMP_INIT1(x1, ctx)
           GR_TMP_INIT2(x1, x2, ctx)
           GR_TMP_INIT3(x1, x2, x3, ctx)
           GR_TMP_INIT4(x1, x2, x3, x4, ctx)
           GR_TMP_INIT5(x1, x2, x3, x4, x5, ctx)

    Allocates one or several temporary elements, all
    initialized to the value 0, assigning the elements
    to the pointers *x1*, *x2*, etc.

.. macro:: GR_TMP_CLEAR1(x1, ctx)
           GR_TMP_CLEAR2(x1, x2, ctx)
           GR_TMP_CLEAR3(x1, x2, x3, ctx)
           GR_TMP_CLEAR4(x1, x2, x3, x4, ctx)
           GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx)

    Corresponding macros to clear temporary variables.

Basic functions
................................................................................

.. function:: int gr_randtest(gr_ptr res, flint_rand_t state, const void * options, gr_ctx_t ctx)

    Sets *res* to a random element of the ring.

.. function:: int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
              int gr_print(gr_srcptr x, gr_ctx_t ctx)
              int gr_println(gr_srcptr x, gr_ctx_t ctx)
              int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx)

    Writes a description of the element *x* to the stream *out*,
    or prints it to *stdout*, or sets *s* to a pointer to
    a heap-allocated string of the description (the user must free
    the string with ``flint_free``). The *println* version prints a
    trailing newline.

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

.. function:: truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* is equal to the element 0, 1 or -1 of the
    ring, respectively.

.. function:: truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Returns whether the elements *x* and *y* are equal.

Arithmetic
........................................................................

User-defined rings should supply ``neg``, ``add``, ``sub``
and ``mul`` methods; the variants with other operand types
have generic fallbacks that may be overridden for performance.
The ``fmpq`` versions may return ``GR_DOMAIN`` if the denominator
is not invertible.

.. function:: int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `-x`.

.. function:: int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to `x + y`.

.. function:: int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to `x - y`.

.. function:: int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to `x \cdot y`.

.. function:: int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to `\mathrm{res } + x \cdot y` or
    `\mathrm{res } - x \cdot y`. Rings may override the default
    implementation to perform this operation in one step without
    allocating a temporary variable.

.. function:: int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `2x`. The default implementation adds *x*
    to itself.

.. function:: int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `x ^ 2`. The default implementation multiplies *x*
    with itself.

Iterated arithmetic operations are best performed using vector
functions.
See in particular :func:`_gr_vec_dot` and :func:`_gr_vec_dot_rev`.

Division
........................................................................

The default implementations of the following methods check for divisors
0, 1, -1 and otherwise return ``GR_UNABLE``.
Particular rings should override the methods when an inversion
or division algorithm is available.
The base rings corresponding to
the following types have complete algorithms
to detect inverses and compute quotients: ``fmpz``, ``fmpq``, ``qqbar``, ``nmod8``.

.. function:: int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to the quotient `x / y` if such an element exists
    in the present ring. Returns the flag ``GR_DOMAIN`` if no such
    quotient exists.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the computation.

    When the ring is not a field, the definition of division may
    vary depending on the ring. A ring implementation may define
    `x / y = x y^{-1}` and return ``GR_DOMAIN`` when `y^{-1}` does not
    exist; alternatively, it may attempt to solve the equation
    `q y = x` (which, for example, gives the usual exact
    division in `\mathbb{Z}`).

.. function:: truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* has a multiplicative inverse in the present ring,
    i.e. whether *x* is a unit.

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

    Sets *res* to the power `x ^ y`, the interpretation of which
    depends on the ring when `y \not \in \mathbb{Z}`.
    Returns the flag ``GR_DOMAIN`` if this power cannot be assigned
    a meaningful value in the present ring, or ``GR_UNABLE`` if
    the implementation is unable to perform the computation.

    For subrings of `\mathbb{C}`, it is implied that the principal
    power `x^y = \exp(y \log(x))` is computed for `x \ne 0`.

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
to detect squares and compute square roots: ``fmpz``, ``fmpq``, ``qqbar``.

In subrings of `\mathbb{C}`, it is implied that the principal
square root is computed; in other cases (e.g. in finite fields),
the choice of root is implementation-dependent.

.. function:: truth_t gr_is_square(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* is a perfect square in the present ring.

.. function:: int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a square root of *x* (respectively reciprocal
    square root) in the present ring, if such an element exists.
    Returns the flag ``GR_DOMAIN`` if *x* is not a perfect square
    (also for zero, when computing the reciprocal square root), or
    ``GR_UNABLE`` if the implementation is unable to perform
    the computation.

Vectors
--------------------------------------------------------------------------------

Low-level vector operations
................................................................................

.. macro:: GR_ENTRY(vec, i, size)

    Macro to access the *i*-th entry of a ``gr_ptr`` or ``gr_srcptr``
    vector *vec*, where each element is ``size`` bytes.

.. function:: void _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)

    Initialize *len* elements of *vec* to the value 0.
    The pointer *vec* must already refer to allocated memory.

.. function:: void _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)

    Clears *len* elements of *vec*.
    This frees memory allocated by individual elements, but
    does not free the memory allocated by *vec* itself.

.. function:: void _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)

    Swap the entries of *vec1* and *vec2*.

.. function:: int _gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, void * options, gr_ctx_t ctx)

.. function:: int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)

.. function:: int _gr_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)

.. function:: int _gr_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: int _gr_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: truth_t _gr_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Memory-managed vectors
................................................................................

.. type:: gr_vec_struct

.. type:: gr_vec_t

.. function:: void gr_vec_init(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Initializes *vec* to a vector of length *len* with elements
    in the ring *ctx*. The length must be nonnegative.
    All entries are set to zero.

.. function:: void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)

    Clears the vector *vec*.

.. function:: gr_ptr gr_vec_entry_ptr(gr_vec_t vec, slong i, gr_ctx_t ctx)

    Returns a pointer to the *i*-th element in the vector *vec*,
    indexed from zero. The index must be in bounds.

.. function:: slong gr_vec_length(const gr_vec_t vec, gr_ctx_t ctx)

    Returns the length of the vector *vec*.

.. function:: void gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Allocates space for at least *len* elements in the vector *vec*.
    This does not change the size of the vector.

.. function:: void gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

    Resizes the vector to length *len*, which must be nonnegative.
    The vector will be extended with zeros.

.. function:: int gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)

    Sets *res* to a copy of the vector *src*.

.. function:: int gr_vec_append(gr_vec_t vec, gr_srcptr x, gr_ctx_t ctx)

    Appends the element *x* to the end of vector *vec*.

Polynomials
--------------------------------------------------------------------------------

See :ref:`gr-poly`.

Implementing rings
--------------------------------------------------------------------------------

.. type:: gr_funcptr

    Typedef for a pointer to a function with signature ``int func(void)``,
    used to represent method table entries.

.. type:: gr_method

    Enumeration type for indexing method tables. Enum values named
    ``GR_METHOD_INIT``,  ``GR_METHOD_ADD_UI``, etc.
    correspond to methods ``gr_init``, ``gr_add_ui``, etc.
    The number of methods is given by ``GR_METHOD_TAB_SIZE``,
    which can be used to declare static method tables.

.. type:: gr_static_method_table

    Typedef for an array of length ``GR_METHOD_TAB_SIZE``
    with :type:`gr_funcptr` entries.

.. function:: int gr_not_implemented(void)

    This function does nothing and returns ``GR_UNABLE``. It is used
    as a generic fallback method when no implementation is available.

.. type:: gr_method_tab_input

    Typedef representing a (index, function pointer) pair.

.. function:: void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)

    Initializes the method table *methods*. This first inserts
    default and generic methods in all slots, and then overwrites
    with the specialized methods listed in *tab*.

