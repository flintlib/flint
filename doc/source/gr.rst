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
        gr_ctx_t ZZ:             /* a parent (context object) */
        gr_ptr x;                /* an element */

        gr_ctx_init_fmpz(ZZ);    /* ZZ = ring of integers with fmpz_t elements */
        GR_TMP_INIT(x, ctx)      /* allocate element on the stack */

        status = gr_set_ui(x, 3, ctx);           /* x = 3 */
        status |= gr_pow_ui(x, x, 100, ctx);     /* x = x ^ 100 */
        status |= gr_println(x, ctx);

        GR_TMP_CLEAR(x, ctx)
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


Mathematical domains
-------------------------------------------------------------------------------

Groups
...............................................................................

.. function:: void gr_ctx_init_perm(gr_ctx_t ctx, ulong n)

    Initializes *ctx* to the symmetric group `S_n` representing
    permutations of `[0, 1, \ldots, n - 1]`.
    Elements are currently represented as pointers (the representation
    may change in the future).

.. function:: void gr_ctx_init_psl2z(gr_ctx_t ctx)

    Initializes *ctx* to the modular group `\text{PSL}(2, \mathbb{Z})`
    with elements of type :type:`psl2z_t`.

.. function:: int gr_ctx_init_dirichlet_group(gr_ctx_t ctx, ulong q)

    Initializes *ctx* to the Dirichlet group `G_q`
    with elements of type :type:`dirichlet_char_t`.
    Fails and returns ``GR_DOMAIN`` if *q* is zero.
    Fails and returns ``GR_UNABLE`` if *q* has a prime factor
    larger than `10^{12}`, which is currently unspported
    by the implementation.

Base rings and fields
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

.. function:: void gr_ctx_init_fmpzi(gr_ctx_t ctx)

    Initializes *ctx* to the ring of Gaussian integers
    `\mathbb{Z}[i]` with elements of type :type:`fmpzi_t`.

.. function:: void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`uint8`. We require `1 \le n \le 255`.

.. function:: void gr_ctx_init_nmod(gr_ctx_t ctx, ulong n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`ulong`. We require `n \ne 0`.

.. function:: void gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`fmpz`. The modulus must be positive.

.. function:: void gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime)

    For a ring initialized with :func:`gr_ctx_init_fmpz_mod`,
    indicate whether the modulus is prime. This can speed up
    some computations.

.. function:: void gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
              void gr_ctx_init_fq_nmod(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
              void gr_ctx_init_fq_zech(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)

    Initializes *ctx* to the finite field `\mathbb{F}_q`
    where `q = p^d`. It is assumed (not checked) that *p* is prime.
    The variable name *var* can be ``NULL`` to use a default.

    The corresponding element types are ``fq_t``, ``fq_nmod_t``, ``fq_zech_t``.
    The ``fq_nmod`` context requires `p < 2^{64}` while ``fq_zech``
    requires `q < 2^{64}` (and in practice a much smaller value
    than this).

.. function:: void gr_ctx_init_real_qqbar(gr_ctx_t ctx)
              void gr_ctx_init_complex_qqbar(gr_ctx_t ctx)

    Initializes *ctx* to the field of real or complex algebraic
    numbers with elements of type :type:`qqbar_t`.

.. function:: void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
              void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the field of real or complex
    numbers represented by elements of type :type:`arb_t`
    and  :type:`acb_t`.

.. function:: void gr_ctx_arb_set_prec(gr_ctx_t ctx, slong prec)
              slong gr_ctx_arb_get_prec(gr_ctx_t ctx)

    Sets or retrieves the bit precision of *ctx* which must be
    an Arb context (this is currently not checked).

.. function:: void gr_ctx_init_real_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_ca(gr_ctx_t ctx)
              void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)

    Initializes *ctx* to the field of real, complex, real algebraic
    or complex algebraic numbers represented by elements of type
    :type:`ca_t`.

.. function:: void gr_ctx_ca_set_option(gr_ctx_t ctx, slong option, slong value)
              slong gr_ctx_ca_get_option(gr_ctx_t ctx, slong option)

    Sets or retrieves options of a Calcium context object.

Floating-point arithmetic
...............................................................................

Although domains of floating-point numbers approximate
real and complex fields, they are not rings or fields.
Floating-point arithmetic can be used in many places where a ring
or field is normally assumed, but predicates like "is field"
return false.

* Equality compares equality of floating-point numbers,
  with the special rule that NaN is not equal to itself.
* In general, the following implementations do not currently
  guarantee correct rounding except for atomic arithmetic operations
  (add, sub, mul, div, sqrt) on real floating-point numbers.

.. function:: void gr_ctx_init_real_float_arf(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the floating-point arithmetic with elements
    of type :type:`arf_t` and a default precision of *prec* bits.

.. function:: void gr_ctx_init_complex_float_acf(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the complex floating-point arithmetic with elements
    of type :type:`acf_t` and a default precision of *prec* bits.

Matrices
...............................................................................

.. function:: void gr_ctx_init_matrix_domain(gr_ctx_t ctx, gr_ctx_t base_ring)

    Initializes *ctx* to the domain of all matrices (of any shape)
    over the given *base_ring*.
    Elements have type :type:`gr_mat_struct`.

.. function:: void gr_ctx_init_matrix_space(gr_ctx_t ctx, gr_ctx_t base_ring, slong n, slong m)

    Initializes *ctx* to the space of matrices over *base_ring*
    with *n* rows and *m* columns.
    Elements have type :type:`gr_mat_struct`.

.. function:: void gr_ctx_init_matrix_ring(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

    Initializes *ctx* to the ring of matrices over *base_ring*
    with *n* rows columns.
    Elements have type :type:`gr_mat_struct`.

Polynomial rings
...............................................................................

.. function:: void gr_ctx_init_polynomial(gr_ctx_t ctx, gr_ctx_t base_ring)

    Initializes *ctx* to a ring of densely represented univariate polynomials
    over the given *base_ring*.
    Elements have type :type:`gr_poly_struct`.

.. function:: void gr_ctx_init_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)

    Initializes *ctx* to a ring of sparsely represented multivariate
    polynomials in *nvars* variables over the given *base_ring*,
    with monomial ordering *ord*.
    Elements have type :type:`gr_mpoly_struct`.

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

Domain properties
...............................................................................

.. function:: truth_t gr_ctx_is_finite(gr_ctx_t ctx)
              truth_t gr_ctx_is_multiplicative_group(gr_ctx_t ctx)
              truth_t gr_ctx_is_ring(gr_ctx_t ctx)
              truth_t gr_ctx_is_commutative_ring(gr_ctx_t ctx)
              truth_t gr_ctx_is_integral_domain(gr_ctx_t ctx)
              truth_t gr_ctx_is_unique_factorization_domain(gr_ctx_t ctx)
              truth_t gr_ctx_is_field(gr_ctx_t ctx)
              truth_t gr_ctx_is_algebraically_closed(gr_ctx_t ctx)
              truth_t gr_ctx_is_finite_characteristic(gr_ctx_t ctx)
              truth_t gr_ctx_is_ordered_ring(gr_ctx_t ctx)

    Returns whether the structure satisfies the respective
    mathematical property.
    The result can be ``T_UNKNOWN``.

.. function:: truth_t gr_ctx_is_exact(gr_ctx_t ctx)

    Returns whether the representation of elements is always exact.

.. function:: truth_t gr_ctx_is_canonical(gr_ctx_t ctx)

    Returns whether the representation of elements is always canonical.

Coercions
...............................................................................

.. function:: int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2)

    Returns 1 if coercing elements into *ctx1* is more meaningful,
    and returns -1 otherwise.

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

.. function:: gr_ptr gr_heap_init(gr_ctx_t ctx)

    Return a pointer to a single new heap-allocated element of *ctx*
    set to 0.

.. function:: void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)

    Free the single heap-allocated element *x* of *ctx* which should
    have been created with :func:`gr_heap_init`.

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

Basic functions
................................................................................

.. function:: int gr_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

    Sets *res* to a random element of the ring.
    The distribution is determined by the ring implementation.
    This normally generates elements non-uniformly to trigger
    corner cases in test code with increased probability.

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

.. function:: int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    Sets *res* to the element *x* of the structure *x_ctx* which
    may be different from *ctx*. This returns the ``GR_DOMAIN`` flag
    if *x* is not an element of *ctx* or cannot be converted
    unambiguously to *ctx*.  The ``GR_UNABLE`` flag is returned
    if the conversion is not implemented.

.. function:: int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx)

    Sets *res* to the string description in *x*.

.. function:: int gr_get_si(slong * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_ui(ulong * res, gr_srcptr x, gr_ctx_t ctx)
              int gr_get_fmpz(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the integer *x*. This returns the ``GR_DOMAIN`` flag
    if *x* is not an integer.

.. function:: int gr_get_fmpq(fmpq_t res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the rational number *x*. This returns the ``GR_DOMAIN``
    flag if *x* is not a rational number.

.. function:: int gr_set_d(gr_ptr res, double x, gr_ctx_t ctx)

    Sets *res* to the value of the floating-point number *x*.
    The interpretation of this conversion depends on the ring.

.. function:: int gr_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)

    Returns a floating-point approximation of *x*. The interpretation
    of this conversion depends on the ring.

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
              int gr_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to `\mathrm{res } + x \cdot y`.
    Rings may override the default
    implementation to perform this operation in one step without
    allocating a temporary variable, without intermediate rounding, etc.

.. function:: int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

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

Iterated arithmetic operations are best performed using vector
functions.
See in particular :func:`_gr_vec_dot` and :func:`_gr_vec_dot_rev`.

.. function:: int _gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
              int _gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
              int _gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the value of the integer polynomial *f* evaluated
    at the argument *x*.

.. function:: int gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)

    Sets *res* to value of the multivariate polynomial *f* (with
    corresponding context object *mctx*) evaluated at the vector
    of arguments in *x*.

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

.. function:: int gr_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_divexact_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_divexact_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_divexact_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_divexact_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)

    Sets *res* to the quotient `x / y`, assuming that this quotient
    is exact in the present ring.
    Rings may optimize this operation by not verifying that the
    division is possible. If the division is not actually exact, the
    implementation may set *res* to a nonsense value and still
    return the ``GR_SUCCESS`` flag.

.. function:: truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx)

    Returns whether *x* has a multiplicative inverse in the present ring,
    i.e. whether *x* is a unit.

.. function:: int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the multiplicative inverse of *x* in the present ring,
    if such an element exists.
    Returns the flag ``GR_DOMAIN`` if *x* is not invertible, or
    ``GR_UNABLE`` if the implementation is unable to perform
    the computation.

Divisibility
........................................................................

.. function:: truth_t gr_divides(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Returns whether *x* divides *y*.

.. function:: int gr_euclidean_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_euclidean_rem(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_euclidean_divrem(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    In a Euclidean ring, these functions perform some version of Euclidean
    division with remainder, where the choice of quotient is
    implementation-defined. For example, it is standard to use
    the round-to-floor quotient in `\mathbb{Z}` and a round-to-nearest quotient in `\mathbb{Z}[i]`.
    In non-Euclidean rings, these functions may implement some generalization of
    Euclidean division with remainder.

.. function:: int gr_gcd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to the greatest common divisor (GCD) of *x* and *y*.

    Since the GCD is only defined uniquely up to multiplication by a unit,
    an implementation-defined representative is chosen.

.. function:: int gr_lcm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to the least common multiple (LCM) of *x* and *y*.

    Since the LCM is only defined uniquely up to multiplication by a unit,
    an implementation-defined representative is chosen.

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

    These methods may return the flag ``GR_DOMAIN`` (or ``GR_UNABLE``)
    when the ring is not a subring of the real or complex numbers.

Ordering methods
........................................................................

.. function:: int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to -1, 0 or 1 according to whether *x* is less than,
    equal or greater than the absolute value of *y*.
    This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

.. function:: int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

    Sets *res* to -1, 0 or 1 according to whether the absolute value
    of *x* is less than, equal or greater than the absolute value of *y*.
    This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

Transcendental functions
........................................................................

.. function:: int gr_pi(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to the constant `\pi`.

.. function:: int gr_euler(gr_ptr res, gr_ctx_t ctx)

    Sets *res* to Euler's constant `\gamma`.

.. function:: int gr_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
              int gr_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_cot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_csc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_exp_pi_i(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_cos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sin_cos_pi(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
              int gr_tan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_cot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_csc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_sinc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sinc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_sinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_cosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sinh_cosh(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)
              int gr_tanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_coth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_sech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_csch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_atan2(gr_ptr res, gr_srcptr y, gr_srcptr x, gr_ctx_t ctx)
              int gr_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_asin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_atan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_asec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acsc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_erf(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_erfc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_erfi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_gamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_lgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_rgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_digamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_zeta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Finite field methods
........................................................................

.. function:: int gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)

.. function:: int gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)

.. function:: int gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)

.. function:: int gr_fq_gen(gr_ptr res, gr_ctx_t ctx)

.. function:: int gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)

.. function:: int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

.. function:: truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

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

.. function:: int _gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)

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

.. function:: int _gr_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)

.. function:: int _gr_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)

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

Required methods
................................................................................

A context object must at minimum define the following methods for a ring:

* init
* clear
* swap
* randtest
* write
* zero
* one
* equal
* set
* set_si
* set_ui
* set_fmpz
* neg
* add
* sub
* mul


Testing rings
--------------------------------------------------------------------------------

.. function:: void gr_test_ring(gr_ctx_t R, slong iters, int test_flags)

    Test correctness of the ring *R*. This calls test functions for
    various methods, each being repeated up to *iters* times.



.. raw:: latex

    \newpage
