.. _fq_default_default:

**fq_default_default.h** -- unified finite fields
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_default_default_ctx_t

    Description.

.. type:: fq_default_default_t

    Description.

Context Management
--------------------------------------------------------------------------------


.. function:: void fq_default_ctx_init(fq_default_ctx_t ctx, const fmpz_t p, slong d, const char *var)

    Initialises the context for prime `p` and extension degree `d`,
    with name ``var`` for the generator.  By default, it will try
    use a Conway polynomial; if one is not available, a random
    irreducible polynomial will be used.

    Assumes that `p` is a prime.

    Assumes that the string ``var`` is a null-terminated string
    of length at least one.

.. function:: void fq_default_ctx_clear(fq_default_ctx_t ctx)

    Clears all memory that has been allocated as part of the context.

.. function:: const fmpz_mod_poly_struct* fq_default_ctx_modulus(const fq_default_ctx_t ctx)

    Returns a pointer to the modulus in the context.

.. function:: slong fq_default_ctx_degree(const fq_default_ctx_t ctx)

    Returns the degree of the field extension
    `[\mathbf{F}_{q} : \mathbf{F}_{p}]`, which
    is equal to `\log_{p} q`.

.. function:: fmpz * fq_default_ctx_prime(const fq_default_ctx_t ctx)

    Returns a pointer to the prime `p` in the context.

.. function:: void fq_default_ctx_order(fmpz_t f, const fq_default_ctx_t ctx)

     Sets `f` to be the size of the finite field.

.. function:: int fq_default_ctx_fprint(FILE * file, const fq_default_ctx_t ctx)

    Prints the context information to ``file``. Returns 1 for a
    success and a negative number for an error.

.. function:: void fq_default_ctx_print(const fq_default_ctx_t ctx)

    Prints the context information to ``stdout``.

.. function:: void fq_default_ctx_randtest(fq_default_ctx_t ctx)

    Initializes ``ctx`` to a random finite field.  Assumes that
    ``fq_default_ctx_init`` has not been called on ``ctx`` already.


Memory management
--------------------------------------------------------------------------------


.. function:: void fq_default_init(fq_default_t rop, const fq_default_ctx_t ctx)

    Initialises the element ``rop``, setting its value to `0`.

.. function:: void fq_default_init2(fq_default_t rop, const fq_default_ctx_t ctx)

    Initialises ``poly`` with at least enough space for it to be an element
    of ``ctx`` and sets it to `0`.

.. function:: void fq_default_clear(fq_default_t rop, const fq_default_ctx_t ctx)

    Clears the element ``rop``.


Basic arithmetic
--------------------------------------------------------------------------------


.. function:: void fq_default_add(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

    Sets ``rop`` to the sum of ``op1`` and ``op2``.

.. function:: void fq_default_sub(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

    Sets ``rop`` to the difference of ``op1`` and ``op2``.

.. function:: void fq_default_sub_one(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

    Sets ``rop`` to the difference of ``op1`` and `1`.

.. function:: void fq_default_neg(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Sets ``rop`` to the negative of ``op``.

.. function:: void fq_default_mul(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``,
    reducing the output in the given context.

.. function:: void fq_default_mul_fmpz(fq_default_t rop, const fq_default_t op, const fmpz_t x, const fq_default_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_default_mul_si(fq_default_t rop, const fq_default_t op, slong x, const fq_default_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_default_mul_ui(fq_default_t rop, const fq_default_t op, ulong x, const fq_default_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_default_sqr(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Sets ``rop`` to the square of ``op``,
    reducing the output in the given context.

.. function:: void fq_default_div(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

    Sets ``rop`` to the quotient of ``op1`` and ``op2``,
    reducing the output in the given context.

.. function:: void fq_default_inv(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Sets ``rop`` to the inverse of the non-zero element ``op``.

.. function:: void fq_default_pow(fq_default_t rop, const fq_default_t op, const fmpz_t e, const fq_default_ctx_t ctx)

    Sets ``rop`` the ``op`` raised to the power `e`.

    Currently assumes that `e \geq 0`.

    Note that for any input ``op``, ``rop`` is set to `1`
    whenever `e = 0`.

.. function:: void fq_default_pow_ui(fq_default_t rop, const fq_default_t op, const ulong e, const fq_default_ctx_t ctx)

    Sets ``rop`` the ``op`` raised to the power `e`.

    Currently assumes that `e \geq 0`.

    Note that for any input ``op``, ``rop`` is set to `1`
    whenever `e = 0`.



Roots
--------------------------------------------------------------------------------


.. function:: int fq_default_sqrt(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

    Sets ``rop`` to the square root of ``op1`` if it is a square, and return
    `1`, otherwise return `0`.

.. function:: void fq_default_pth_root(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

    Sets ``rop`` to a `p^{th}` root root of ``op1``.  Currently,
    this computes the root by raising ``op1`` to `p^{d-1}` where
    `d` is the degree of the extension.

.. function:: int fq_default_is_square(const fq_default_t op, const fq_default_ctx_t ctx)

    Return ``1`` if ``op`` is a square.

Output
--------------------------------------------------------------------------------


.. function:: int fq_default_fprint_pretty(FILE *file, const fq_default_t op, const fq_default_ctx_t ctx)

    Prints a pretty representation of ``op`` to ``file``.

    In the current implementation, always returns `1`.  The return code is
    part of the function's signature to allow for a later implementation to
    return the number of characters printed or a non-positive error code.

.. function:: int fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx)

    Prints a pretty representation of ``op`` to ``stdout``.

    In the current implementation, always returns `1`.  The return code is
    part of the function's signature to allow for a later implementation to
    return the number of characters printed or a non-positive error code.

.. function:: void fq_default_fprint(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)

    Prints a representation of ``op`` to ``file``.

.. function:: void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx)

    Prints a representation of ``op`` to ``stdout``.

.. function:: char * fq_default_get_str(const fq_default_t op, const fq_default_ctx_t ctx)

    Returns the plain FLINT string representation of the element
    ``op``.

.. function:: char * fq_default_get_str_pretty(const fq_default_t op, const fq_default_ctx_t ctx)

    Returns a pretty representation of the element ``op`` using the
    null-terminated string ``x`` as the variable name.


Randomisation
--------------------------------------------------------------------------------


.. function:: void fq_default_randtest(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

    Generates a random element of `\mathbf{F}_q`.

.. function:: void fq_default_randtest_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

    Generates a random non-zero element of `\mathbf{F}_q`.

.. function:: void fq_default_rand(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

    Generates a high quality random element of `\mathbf{F}_q`.

.. function:: void fq_default_rand_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

    Generates a high quality non-zero random element of `\mathbf{F}_q`.


Assignments and conversions
--------------------------------------------------------------------------------


.. function:: void fq_default_set(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Sets ``rop`` to ``op``.

.. function:: void fq_default_set_si(fq_default_t rop, const slong x, const fq_default_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_default_set_ui(fq_default_t rop, const ulong x, const fq_default_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_default_set_fmpz(fq_default_t rop, const fmpz_t x, const fq_default_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_default_swap(fq_default_t op1, fq_default_t op2, const fq_default_ctx_t ctx)

    Swaps the two elements ``op1`` and ``op2``.

.. function:: void fq_default_zero(fq_default_t rop, const fq_default_ctx_t ctx)

    Sets ``rop`` to zero.

.. function:: void fq_default_one(fq_default_t rop, const fq_default_ctx_t ctx)

    Sets ``rop`` to one, reduced in the given context.

.. function:: void fq_default_gen(fq_default_t rop, const fq_default_ctx_t ctx)

    Sets ``rop`` to a generator for the finite field.
    There is no guarantee this is a multiplicative generator of
    the finite field.


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_default_is_zero(const fq_default_t op, const fq_default_ctx_t ctx)

    Returns whether ``op`` is equal to zero.

.. function:: int fq_default_is_one(const fq_default_t op, const fq_default_ctx_t ctx)

    Returns whether ``op`` is equal to one.

.. function:: int fq_default_equal(const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

    Returns whether ``op1`` and ``op2`` are equal.


Special functions
--------------------------------------------------------------------------------


.. function:: void fq_default_trace(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Sets ``rop`` to the trace of ``op``.

    For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the
    trace of `a` as the trace of this map.  Equivalently, if `\Sigma`
    generates `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of
    `a` is equal to `\sum_{i=0}^{d-1} \Sigma^i (a)`, where `d =
    \log_{p} q`.

.. function:: void fq_default_norm(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

    Computes the norm of ``op``.

    For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the norm
    of `a` as the determinant of this map.  Equivalently, if `\Sigma` generates
    `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of `a` is equal to
    `\prod_{i=0}^{d-1} \Sigma^i (a)`, where
    `d = \text{dim}_{\mathbf{F}_p}(\mathbf{F}_q)`.

    Algorithm selection is automatic depending on the input.

.. function:: void fq_default_frobenius(fq_default_t rop, const fq_default_t op, slong e, const fq_default_ctx_t ctx)

    Evaluates the homomorphism `\Sigma^e` at ``op``.

    Recall that `\mathbf{F}_q / \mathbf{F}_p` is Galois with Galois group
    `\langle \sigma \rangle`, which is also isomorphic to
    `\mathbf{Z}/d\mathbf{Z}`, where
    `\sigma \in \operatorname{Gal}(\mathbf{F}_q/\mathbf{F}_p)` is the Frobenius element
    `\sigma \colon x \mapsto x^p`.

