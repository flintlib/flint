.. _gr:

**gr.h** -- generic structures and their elements
===============================================================================

Introduction
-------------------------------------------------------------------------------

Parents and elements
...............................................................................

To work with an element `x \in R` of a particular mathematical
structure *R*, we use a context object to represent *R*
(the "parent" of `x`). Elements are passed around as pointers.
Note:

* Parents are not stored as part of the elements; the user must
  track the context objects for all variables.
* Operations are strictly type-stable:
  elements only change parent when performing an explicit conversion.

The structure *R* will typically be a *ring*, but the framework supports
general
objects (including groups, monoids, and sets without any particular
structure whatsoever). We use these terms in a strict mathematical
sense: a "ring" must exactly satisfy the ring axioms.
It can have inexact *representations*, but this inexactness
must be handled rigorously.

To give an idea of how the interface works, this example program
computes `3^{100}` in the ring of integers and prints the value::

    #include "gr.h"

    int main()
    {
        int status;
        gr_ctx_t ZZ;             /* a parent (context object) */
        gr_ptr x;                /* an element */

        gr_ctx_init_fmpz(ZZ);    /* ZZ = ring of integers with fmpz_t elements */
        GR_TMP_INIT(x, ZZ)       /* allocate element on the stack */

        status = gr_set_ui(x, 3, ZZ);           /* x = 3 */
        status |= gr_pow_ui(x, x, 100, ZZ);     /* x = x ^ 100 */
        status |= gr_println(x, ZZ);

        GR_TMP_CLEAR(x, ZZ)
        gr_ctx_clear(ZZ);

        return status;
    }

Parent and element types
...............................................................................

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

    A context object representing a mathematical structure *R*.
    It contains the following data:

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

There is no type to represent a single generic element
as a struct since we do not know the size of a generic element at
compile time.
Memory for single elements can either be allocated on the stack
with the special macros provided below, or as usual with ``malloc``.
Methods can also be used with particular C types like ``fmpz_t``
when the user knows the type.
Users may wish to define their own union types when only some
particular types will appear in an application.

Error handling
...............................................................................

To compute over a structure `R`, it is useful to conceptually extend
to a larger set `R' = R \cup \{ \text{undefined}, \text{unknown} \}`.

* Adding an *undefined* (error) value allows us to extend partial functions
  to total functions.
* An *unknown* value is useful in cases where a result
  may exist in principle but cannot be computed.

An alternative to having an *undefined* value
is to choose some arbitrary default value in `R`,
say `\text{undefined} = 0` in a ring. This is often done in
proof assistants, but in a regular programming environment,
we typically want some way to detect domain errors.

Representing `R'` as a type-level extension of `R` is tricky in C
since we would either have to
wrap elements in a larger structure
or reserve bit patterns in each type for special values.
In any case, it is useful to assume in low-level code
that elements *really represent elements of the intended structure*
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

.. macro:: GR_TEST_FAIL

    Test failure. This is only used in test code.

When the status code is any other value than ``GR_SUCCESS``, any
output variables may be set to meaningless values.

C functions that return a status code are marked with the
``WARN_UNUSED_RESULT`` attribute. This allows compilers to
emit warnings when the status code is ignored.

Flags can be OR'ed and checked only at the top level of a computation
to avoid complex control flow::

    status = GR_SUCCESS;
    gr |= gr_add(res, a, b, ctx);
    gr |= gr_pow_ui(res, res, 2, ctx);
    ...

If we do not care about recovering from *undefined*/*unknown* results,
the following macro is useful:

.. macro:: GR_MUST_SUCCEED(expr)

    Evaluates *expr* and asserts that the return value is
    ``GR_SUCCESS``. On failure, calls ``flint_abort()``.

For uniformity, most operations return a status code, even operations
that are not typically expected to fail. Exceptions include the
following:

* Pure "container" operations like ``init``, ``clear`` and ``swap``
  do not return a status code.

* Pure predicate functions (see below)
  return ``T_TRUE`` / ``T_FALSE`` / ``T_UNKNOWN``
  instead of computing a separate boolean value and error code.


Predicates
...............................................................................

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


Context operations
-------------------------------------------------------------------------------

.. function:: slong gr_ctx_sizeof_elem(gr_ctx_t ctx)

    Return ``sizeof(type)`` where ``type`` is the underlying C
    type for elements of *ctx*.

.. function:: int gr_ctx_clear(gr_ctx_t ctx)

    Clears the context object *ctx*, freeing any memory
    allocated by this object.

    Some context objects may require that no elements are cleared after calling
    this method, and may leak memory if not all elements have
    been cleared when calling this method.

    If *ctx* is derived from a base ring, the base ring context
    may also be required to stay alive until after this
    method is called.

.. function:: int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx)
              int gr_ctx_print(gr_ctx_t ctx)
              int gr_ctx_println(gr_ctx_t ctx)
              int gr_ctx_get_str(char ** s, gr_ctx_t ctx)

    Writes a description of the structure *ctx* to the stream *out*,
    prints it to *stdout*, or sets *s* to a pointer to
    a heap-allocated string of the description (the user must free
    the string with ``flint_free``).
    The *println* version prints a trailing newline.

.. function:: int gr_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
              int gr_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)

    Set the name of the generator (univariate polynomial ring,
    finite field, etc.) or generators (multivariate).
    The name is used when printing and may be used to choose
    coercions.

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

.. function:: void gr_set_shallow(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a shallow copy of *x*, copying the struct data.

.. function:: gr_ptr gr_heap_init(gr_ctx_t ctx)

    Return a pointer to a single new heap-allocated element of *ctx*
    set to 0.

.. function:: void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)

    Free the single heap-allocated element *x* of *ctx* which should
    have been created with :func:`gr_heap_init`.

.. function:: gr_ptr gr_heap_init_vec(slong len, gr_ctx_t ctx)

    Return a pointer to a new heap-allocated vector of *len*
    initialized elements.

.. function:: void gr_heap_clear_vec(gr_ptr x, slong len, gr_ctx_t ctx)

    Clear the *len* elements in the heap-allocated vector *len* and
    free the vector itself.

The following macros support allocating temporary variables efficiently.
Data will be allocated on the stack using ``alloca`` unless
the size is excessive (risking stack overflow), in which case
the implementation transparently switches to ``malloc``/``free``
instead. The usage pattern is as follows::

    {
        gr_ptr x, y;
        GR_TMP_INIT2(x1, x2, ctx);

        /* do computations with x1, x2 */

        GR_TMP_CLEAR2(x1, x2, ctx);
    }

Init and clear macros must match exactly, as variables may be
allocated contiguously in a block.

*Warning:* never use these macros directly inside a loop.
This is likely to overflow the stack, as memory will not
be reclaimed until the function exits.
Instead, allocate the needed space before entering
any loops, move the loop body to a separate function,
or allocate the memory on the heap if needed.

.. macro:: GR_TMP_INIT_VEC(vec, len, ctx)
           GR_TMP_CLEAR_VEC(vec, len, ctx)

    Allocates and frees a vector of *len* contiguous elements, all
    initialized to the value 0, assigning the first element
    to the pointer *vec*.

.. macro:: GR_TMP_INIT(x1, ctx)
           GR_TMP_INIT2(x1, x2, ctx)
           GR_TMP_INIT3(x1, x2, x3, ctx)
           GR_TMP_INIT4(x1, x2, x3, x4, ctx)
           GR_TMP_INIT5(x1, x2, x3, x4, x5, ctx)

    Allocates one or several temporary elements, all
    initialized to the value 0, assigning the elements
    to the pointers *x1*, *x2*, etc.

.. macro:: GR_TMP_CLEAR(x1, ctx)
           GR_TMP_CLEAR2(x1, x2, ctx)
           GR_TMP_CLEAR3(x1, x2, x3, ctx)
           GR_TMP_CLEAR4(x1, x2, x3, x4, ctx)
           GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx)

    Corresponding macros to clear temporary variables.

Random elements
................................................................................

.. function:: int gr_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

    Sets *res* to a random element of the domain *ctx*.
    The distribution is determined by the implementation.
    Typically the distribution is non-uniform in order to
    find corner cases more easily in test code.

.. function:: int gr_randtest_not_zero(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

    Sets *res* to a random nonzero element of the domain *ctx*.
    This operation will fail and return ``GR_DOMAIN`` in the zero ring.

.. function:: int gr_randtest_small(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

    Sets *res* to a "small" element of the domain *ctx*.
    This is suitable for randomized testing where a "large" argument
    could result in excessive computation time.

Input, output and string conversion
................................................................................

.. function:: int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
              int gr_print(gr_srcptr x, gr_ctx_t ctx)
              int gr_println(gr_srcptr x, gr_ctx_t ctx)
              int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx)

    Writes a description of the element *x* to the stream *out*,
    or prints it to *stdout*, or sets *s* to a pointer to
    a heap-allocated string of the description (the user must free
    the string with ``flint_free``). The *println* version prints a
    trailing newline.

.. function:: int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx)

    Sets *res* to the string description in *x*.

.. function:: int gr_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)
              int gr_get_str_n(char ** s, gr_srcptr x, slong n, gr_ctx_t ctx)

    String conversion where real and complex numbers may be rounded
    to *n* digits.

Assignment and conversions
................................................................................

.. function:: int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a copy of the element *x*.

.. function:: int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    Sets *res* to the element *x* of the structure *x_ctx* which
    may be different from *ctx*. This returns the ``GR_DOMAIN`` flag
    if *x* is not an element of *ctx* or cannot be converted
    unambiguously to *ctx*.  The ``GR_UNABLE`` flag is returned
    if the conversion is not implemented.

.. function:: int gr_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
              int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx)
              int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
              int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx)
              int gr_set_d(gr_ptr res, double x, gr_ctx_t ctx)

    Sets *res* to the value *x*. If no reasonable conversion to the
    domain *ctx* is possible, returns ``GR_DOMAIN``.

.. function:: int gr_get_si(slong * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_ui(ulong * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_fmpz(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_fmpq(fmpq_t res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the value *x*. This returns the ``GR_DOMAIN`` flag
    if *x* cannot be converted to the target type.
    For floating-point output types, the output may be rounded.

.. function:: int gr_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t a, const fmpz_t b, gr_ctx_t ctx)
              int gr_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_srcptr x, gr_ctx_t ctx)

    Set or retrieve a dyadic number `a \cdot 2^b`.

.. function:: int gr_set_fmpz_10exp_fmpz(gr_ptr res, const fmpz_t a, const fmpz_t b, gr_ctx_t ctx)

    Set to a decimal number `a \cdot 10^b`.

.. function:: int gr_get_fexpr(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_fexpr_serialize(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to a symbolic expression representing *x*.
    The *serialize* version may generate a representation of the
    internal representation which is not intended to be human-readable.

.. function:: int gr_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t x, gr_ctx_t ctx)

    Sets *res* to the evaluation of the expression *x* in the
    given ring or structure.
    The user must provide vectors *inputs* and *outputs* which
    may be empty initially and which may be used as scratch space
    during evaluation. Non-empty vectors may be given to map symbols
    to predefined values.

Special values
................................................................................

.. function:: int gr_zero(gr_ptr res, gr_ctx_t ctx)
              int gr_one(gr_ptr res, gr_ctx_t ctx)
              int gr_neg_one(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to the ring element 0, 1 or -1.

.. function:: int gr_gen(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to a generator of this domain. The meaning of
    "generator" depends on the domain.

.. function:: int gr_gens(gr_vec_t res, gr_ctx_t ctx)
              int gr_gens_recursive(gr_vec_t res, gr_ctx_t ctx)

    Sets *res* to a vector containing the generators of this domain
    where this makes sense, for example in a multivariate polynomial
    ring. The *recursive* version also includes any generators
    of the base ring, and of any recursive base rings.

Basic properties
........................................................................

.. function:: truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* is equal to the ring element 0, 1 or -1,
    respectively.

.. function:: truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Returns whether the elements *x* and *y* are equal.

.. function:: truth_t gr_is_integer(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* represents an integer.

.. function:: truth_t gr_is_rational(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* represents a rational number.

Arithmetic
........................................................................

User-defined rings should supply ``neg``, ``add``, ``sub``
and ``mul`` methods; the variants with other operand types
have generic fallbacks that may be overridden for performance.
The ``fmpq`` versions may return ``GR_DOMAIN`` if the denominator
is not invertible.
The *other* versions accept operands belonging to a different domain,
attempting to perform a coercion into the target domain.

.. function:: int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `-x`.

.. function:: int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to `x + y`.

.. function:: int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to `x - y`.

.. function:: int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to `x \cdot y`.

.. function:: int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

    Sets *res* to `\mathrm{res } + x \cdot y`.
    Rings may override the default
    implementation to perform this operation in one step without
    allocating a temporary variable, without intermediate rounding, etc.

.. function:: int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

    Sets *res* to `\mathrm{res } - x \cdot y`.
    Rings may override the default
    implementation to perform this operation in one step without
    allocating a temporary variable, without intermediate rounding, etc.

.. function:: int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `2x`. The default implementation adds *x*
    to itself.

.. function:: int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to `x ^ 2`. The default implementation multiplies *x*
    with itself.

.. function:: int gr_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)

    Sets *res* to `x \cdot 2^y`. This may perform `x \cdot 2^{-y}`
    when *y* is negative, allowing exact division by powers of two
    even if `2^{y}` is not representable.

Iterated arithmetic operations are best performed using vector
functions.
See in particular :func:`_gr_vec_dot` and :func:`_gr_vec_dot_rev`.

Division
........................................................................

The default implementations of the following methods check for divisors
0, 1, -1 and otherwise return ``GR_UNABLE``.
Particular rings should override the methods when an inversion
or division algorithm is available.

.. function:: truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* has a multiplicative inverse in the present ring,
    i.e. whether *x* is a unit.

.. function:: int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the multiplicative inverse of *x* in the present ring,
    if such an element exists.
    Returns the flag ``GR_DOMAIN`` if *x* is not invertible, or
    ``GR_UNABLE`` if the implementation is unable to perform
    the computation.

.. function:: int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to the quotient `x / y`. In a field, this returns
    ``GR_DOMAIN`` if `y` is zero; in an integral domain,
    it returns ``GR_DOMAIN`` if `y` is zero or if the quotient
    does not exist. In a non-integral domain, we consider a quotient
    to exist only if it is unique, and otherwise return ``GR_DOMAIN``;
    see :func:`gr_div_nonunique` for a different behavior.

    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the computation.

.. function:: int gr_div_nonunique(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to an arbitrary solution `q` of the equation `x = q y`.
    Returns the flag ``GR_DOMAIN`` if no such solution exists.
    Returns the flag ``GR_UNABLE`` if the implementation is unable
    to perform the computation.
    This method allows dividing `x / y` in some cases where :func:`gr_div` fails:

    * `0 / 0` has solutions (for example, 0) in any ring.
    * It allows solving division problems in nonintegral domains.
      For example, it allows assigning a value to `6 / 2` in
      `R = \mathbb{Z}/10\mathbb{Z}` even though `2^{-1}` does not exist
      in `R`. In this case, both 3 and 8 are possible solutions, and which
      one is chosen is unpredictable.

.. function:: truth_t gr_divides(gr_srcptr d, gr_srcptr x, gr_ctx_t ctx)

    Returns whether `d \mid x`; that is, whether there is an element `q`
    such that `x = dq`. Note that this corresponds to divisibility
    in the sense of :func:`gr_div_nonunique`, which is weaker than that
    of :func:`gr_div`. For example, `0 \mid 0` is true even
    in rings where `0 / 0` is undefined.

.. function:: int gr_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_divexact_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_divexact_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_divexact_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_divexact_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_divexact(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to the quotient `x / y`, assuming that this quotient
    is exact in the present ring.
    Rings may optimize this operation by not verifying that the
    division is possible. If the division is not actually exact, the
    implementation may set *res* to a nonsense value and still
    return the ``GR_SUCCESS`` flag.

.. function:: int gr_euclidean_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_euclidean_rem(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_euclidean_divrem(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    In a Euclidean ring, these functions perform some version of Euclidean
    division with remainder, where the choice of quotient is
    implementation-defined. For example, it is standard to use
    the round-to-floor quotient in `\mathbb{Z}` and a round-to-nearest quotient in `\mathbb{Z}[i]`.
    In non-Euclidean rings, these functions may implement some generalization of
    Euclidean division with remainder.

Powering
........................................................................

.. function:: int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

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

Greatest common divisors
........................................................................

.. function:: int gr_gcd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to a greatest common divisor (GCD) of *x* and *y*.
    Since the GCD is unique only up to multiplication by a unit,
    an implementation-defined representative is chosen.

.. function:: int gr_lcm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to a least common multiple (LCM) of *x* and *y*.
    Since the LCM is unique only up to multiplication by a unit,
    an implementation-defined representative is chosen.

Factorization
........................................................................

.. function:: int gr_factor(gr_ptr c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)

    Given `x \in R`, computes a factorization

        `x = c {f_1}^{e_1} \ldots {f_n}^{e_n}`

    where `f_k` will be irreducible or prime (depending on `R`).

    The prefactor `c` stores a unit, sign, or coefficient, e.g.\ the
    sign `-1`, `0` or `+1` in `\mathbb{Z}`, or a sign multiplied
    by the coefficient content in `\mathbb{Z}[x]`.
    Note that this function outputs `c` as an element of the
    same ring as the input: for example, in `\mathbb{Z}[x]`,
    `c` will be a constant polynomial rather than an
    element of the coefficient ring.
    The exponents `e_k` are output as a vector of ``fmpz`` elements.

    The factors `f_k` are guaranteed to be distinct,
    but they are not guaranteed to be sorted in any particular
    order.

Fractions
........................................................................

.. function:: int gr_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Return a numerator `p` and denominator `q` such that `x = p/q`.
    For typical fraction fields, the denominator will be minimal
    and canonical.
    However, some rings may return an arbitrary denominator as long
    as the numerator matches.
    The default implementations simply return `p = x` and `q = 1`.

Integer and complex parts
........................................................................

.. function:: int gr_floor(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_ceil(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_trunc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_nint(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    In the real and complex numbers, sets *res* to the integer closest
    to *x*, respectively rounding towards minus infinity, plus infinity,
    zero, or the nearest integer (with tie-to-even).

.. function:: int gr_abs(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the absolute value of *x*, which maybe defined
    both in complex rings and in any ordered ring.

.. function:: int gr_i(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to the imaginary unit.

.. function:: int gr_conj(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_re(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_im(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_csgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_arg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    These methods may return the flag ``GR_DOMAIN`` (or ``GR_UNABLE``)
    when the ring is not a subring of the real or complex numbers.

Infinities and extended values
........................................................................

.. function:: int gr_pos_inf(gr_ptr res, gr_ctx_t ctx)
              int gr_neg_inf(gr_ptr res, gr_ctx_t ctx)
              int gr_uinf(gr_ptr res, gr_ctx_t ctx)
              int gr_undefined(gr_ptr res, gr_ctx_t ctx)
              int gr_unknown(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to the signed positive infinity `+\infty`, signed negative infinity `-\infty`, unsigned infinity `{\tilde \infty}`, *Undefined*, or *Unknown*, respectively.

Ordering methods
........................................................................

.. function:: int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

    Sets *res* to -1, 0 or 1 according to whether *x* is less than,
    equal or greater than *y*.
    This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

.. function:: int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

    Sets *res* to -1, 0 or 1 according to whether the absolute value
    of *x* is less than, equal or greater than the absolute value of *y*.
    This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

.. function:: truth_t gr_le(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_lt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_ge(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_gt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_abs_le(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_abs_lt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_abs_ge(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              truth_t gr_abs_gt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Wrappers of ``gr_cmp`` and ``gr_cmpabs`` returning truth values
    for the comparison operations ``<=``, ``<``, ``>=``, ``>``.

.. function:: int gr_min(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_max(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Minimum and maximum value.

Enclosure and interval methods
........................................................................

.. function:: int gr_set_interval_mid_rad(gr_ptr res, gr_srcptr m, gr_srcptr r, gr_ctx_t ctx)

    In ball representations of the real numbers, sets *res* to
    the interval `m \pm r`.

    In vector spaces over the real numbers represented using balls,
    intervals are handled independently for the generators;
    for example, in the complex numbers, `a + b i \pm (0.1 + 0.2 i)`
    is equivalent to `(a \pm 0.1) + (b \pm 0.2) i`.

Finite field methods
........................................................................

.. function:: int gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)

.. function:: int gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)

.. function:: int gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)

.. function:: int gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)

.. function:: int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)



.. raw:: latex

    \newpage
