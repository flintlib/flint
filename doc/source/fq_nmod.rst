.. _fq-nmod:

**fq_nmod.h** -- finite fields (word-size characteristic)
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_nmod_ctx_struct

.. type:: fq_nmod_ctx_t

    Description.

.. type:: fq_nmod_struct

.. type:: fq_nmod_t

    Description.

Context Management
--------------------------------------------------------------------------------


.. function:: void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)

    Initialises the context for prime `p` and extension degree `d`,
    with name ``var`` for the generator.  By default, it will try
    use a Conway polynomial; if one is not available, a random
    irreducible polynomial will be used.

    Assumes that `p` is a prime.

    Assumes that the string ``var`` is a null-terminated string
    of length at least one.

.. function:: int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)

    Attempts to initialise the context for prime `p` and extension
    degree `d`, with name ``var`` for the generator using a Conway
    polynomial for the modulus.

    Returns `1` if the Conway polynomial is in the database for the
    given size and the initialization is successful; otherwise,
    returns `0`.

    Assumes that `p` is a prime.

    Assumes that the string ``var`` is a null-terminated string
    of length at least one.

.. function:: void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)

    Initialises the context for prime `p` and extension degree `d`,
    with name ``var`` for the generator using a Conway polynomial
    for the modulus.

    Assumes that `p` is a prime.

    Assumes that the string ``var`` is a null-terminated string
    of length at least one.

.. function:: void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx, nmod_poly_t modulus, const char *var)

    Initialises the context for given ``modulus`` with name
    ``var`` for the generator.

    Assumes that ``modulus`` is an irreducible polynomial over
    `\mathbf{F}_{p}`.

    Assumes that the string ``var`` is a null-terminated string
    of length at least one.

.. function:: void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx)

    Clears all memory that has been allocated as part of the context.

.. function:: const nmod_poly_struct* fq_nmod_ctx_modulus(const fq_nmod_ctx_t ctx)

    Returns a pointer to the modulus in the context.

.. function:: long fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)

    Returns the degree of the field extension
    `[\mathbf{F}_{q} : \mathbf{F}_{p}]`, which
    is equal to `\log_{p} q`.

.. function:: fmpz * fq_nmod_ctx_prime(const fq_nmod_ctx_t ctx)

    Returns a pointer to the prime `p` in the context.

.. function:: void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)

     Sets `f` to be the size of the finite field.

.. function:: int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx)

    Prints the context information to ``file``. Returns 1 for a
    success and a negative number for an error.

.. function:: void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)

    Prints the context information to {``stdout``.

.. function:: void fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx)

    Initializes ``ctx`` to a random finite field.  Assumes that
    ``fq_nmod_ctx_init`` has not been called on ``ctx`` already.


.. function:: void fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx)

    Initializes ``ctx`` to a random extension of a word-sized prime
    field.  The modulus may or may not be irreducible.  Assumes that
    ``fq_nmod_ctx_init`` as not been called on ``ctx`` already.


Memory management
--------------------------------------------------------------------------------


.. function:: void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Initialises the element ``rop``, setting its value to `0`. Currently, the behaviour is identical to ``fq_nmod_init2``, as it also ensures ``rop`` has enough space for it to be an element of ``ctx``, this may change in the future.

.. function:: void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Initialises ``rop`` with at least enough space for it to be an element
    of ``ctx`` and sets it to `0`.

.. function:: void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Clears the element ``rop``.

.. function:: void _fq_nmod_sparse_reduce(mp_ptr R, slong lenR, const fq_nmod_ctx_t ctx)

    Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    modulus of ``ctx``.

.. function:: void _fq_nmod_dense_reduce(mp_ptr R, slong lenR, const fq_nmod_ctx_t ctx)

    Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    modulus of ``ctx`` using Newton division.

.. function:: void _fq_nmod_reduce(mp_ptr r, slong lenR, const fq_nmod_ctx_t ctx)

    Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    modulus of ``ctx``.  Does either sparse or dense reduction
    based on ``ctx->sparse_modulus``.

.. function:: void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Reduces the polynomial ``rop`` as an element of
    `\mathbf{F}_p[X] / (f(X))`.


Basic arithmetic
--------------------------------------------------------------------------------


.. function:: void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the sum of ``op1`` and ``op2``.

.. function:: void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the difference of ``op1`` and ``op2``.

.. function:: void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the difference of ``op1`` and `1`.

.. function:: void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the negative of ``op``.

.. function:: void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``,
    reducing the output in the given context.

.. function:: void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_nmod_mul_ui(fq_nmod_t rop, const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` and `x`,
    reducing the output in the given context.

.. function:: void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the square of ``op``,
    reducing the output in the given context.

.. function:: void _fq_nmod_inv(mp_ptr *rop, mp_srcptr *op, slong len, const fq_nmod_ctx_t ctx)

    Sets ``(rop, d)`` to the inverse of the non-zero element
    ``(op, len)``.

.. function:: void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the inverse of the non-zero element ``op``.

.. function:: void fq_nmod_gcdinv(fq_nmod_t f, fq_nmod_t inv, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

     Sets ``inv`` to be the inverse of ``op`` modulo the modulus
     of ``ctx``.  If ``op`` is not invertible, then ``f`` is
     set to a factor of the modulus; otherwise, it is set to one.

.. function:: void _fq_nmod_pow(mp_ptr *rop, mp_srcptr *op, slong len, const fmpz_t e, const fq_nmod_ctx_t ctx)

    Sets ``(rop, 2*d-1)`` to ``(op,len)`` raised to the power `e`,
    reduced modulo `f(X)`, the modulus of ``ctx``.

    Assumes that `e \geq 0` and that ``len`` is positive and at most `d`.

    Although we require that ``rop`` provides space for
    `2d - 1` coefficients, the output will be reduced modulo
    `f(X)`, which is a polynomial of degree `d`.

    Does not support aliasing.

.. function:: void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t e, const fq_nmod_ctx_t ctx)

    Sets ``rop`` the ``op`` raised to the power `e`.

    Currently assumes that `e \geq 0`.

    Note that for any input ``op``, ``rop`` is set to `1`
    whenever `e = 0`.

.. function:: void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op, const ulong e, const fq_nmod_ctx_t ctx)

    Sets ``rop`` the ``op`` raised to the power `e`.

    Currently assumes that `e \geq 0`.

    Note that for any input ``op``, ``rop`` is set to `1`
    whenever `e = 0`.


Roots
--------------------------------------------------------------------------------


.. function:: void fq_nmod_sqrt(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)                                                
    Sets ``rop`` to the square root of ``op1`` if it is a square, and return
    `1`, otherwise return `0`.

.. function:: void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to a `p^{th}` root root of ``op1``.  Currently,
    this computes the root by raising ``op1`` to `p^{d-1}` where
    `d` is the degree of the extension.

.. function:: int fq_nmod_is_square(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Return ``1`` if ``op`` is a square.

Output
--------------------------------------------------------------------------------


.. function:: int fq_nmod_fprint_pretty(FILE *file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Prints a pretty representation of ``op`` to ``file``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Prints a pretty representation of ``op`` to ``stdout``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: void fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Prints a representation of ``op`` to ``file``.

    For further details on the representation used, see
    ``nmod_poly_fprint()``.

.. function:: void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Prints a representation of ``op`` to ``stdout``.

    For further details on the representation used, see
    ``nmod_poly_print()``.

.. function:: char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns the plain FLINT string representation of the element
    ``op``.

.. function:: char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns a pretty representation of the element ``op`` using the
    null-terminated string ``x`` as the variable name.


Randomisation
--------------------------------------------------------------------------------


.. function:: void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Generates a random element of `\mathbf{F}_q`.

.. function:: void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Generates a random non-zero element of `\mathbf{F}_q`.

.. function:: void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Generates a random element of `\mathbf{F}_q` which has an
    underlying polynomial with dense coefficients.

.. function:: void fq_nmod_rand(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)

    Generates a high quality random element of `\mathbf{F}_q`.


Assignments and conversions
--------------------------------------------------------------------------------


.. function:: void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to ``op``.

.. function:: void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to ``x``, considered as an element of
    `\mathbf{F}_p`.

.. function:: void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx)

    Swaps the two elements ``op1`` and ``op2``.

.. function:: void fq_nmod_zero(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to zero.

.. function:: void fq_nmod_one(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to one, reduced in the given context.

.. function:: void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to a generator for the finite field.
    There is no guarantee this is a multiplicative generator of
    the finite field.

.. function:: void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx);

    Set ``a`` to a representative of ``b`` in ``ctx``.
    The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

.. function:: void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b, const fq_nmod_ctx_t ctx);

    Set ``a`` to the element in ``ctx`` with representative ``b``.
    The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

.. function:: void fq_nmod_get_nmod_mat(nmod_mat_t col, const fq_nmod_t a, const fq_nmod_ctx_t ctx)

    Convert ``a`` to a column vector of length ``degree(ctx)``.

.. function:: void fq_nmod_set_nmod_mat(fq_nmod_t a, const nmod_mat_t col, const fq_nmod_ctx_t ctx)

    Convert a column vector ``col`` of length ``degree(ctx)`` to an element of ``ctx``.


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns whether ``op`` is equal to zero.

.. function:: int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns whether ``op`` is equal to one.

.. function:: int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)

    Returns whether ``op1`` and ``op2`` are equal.

.. function:: int fq_nmod_is_invertible(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns whether ``op`` is an invertible element.

.. function:: int fq_nmod_is_invertible_f(fq_nmod_t f, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns whether ``op`` is an invertible element.  If it is not,
    then ``f`` is set of a factor of the modulus.

.. function:: int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)

    Return ``1`` (resp. ``-1``, or ``0``) if ``a`` is after (resp. before, same as) ``b`` in some arbitrary but fixed total ordering of the elements.


Special functions
--------------------------------------------------------------------------------


.. function:: void _fq_nmod_trace(fmpz_t rop, mp_srcptr *op, slong len, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the trace of the non-zero element ``(op, len)``
    in `\mathbf{F}_{q}`.

.. function:: void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the trace of ``op``.

    For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the
    trace of `a` as the trace of this map.  Equivalently, if `\Sigma`
    generates `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of
    `a` is equal to `\sum_{i=0}^{d-1} \Sigma^i (a)`, where `d =
    \log_{p} q`.

.. function:: void _fq_nmod_norm(fmpz_t rop, mp_srcptr *op, slong len, const fq_nmod_ctx_t ctx)

    Sets ``rop`` to the norm of the non-zero element ``(op, len)``
    in `\mathbf{F}_{q}`.

.. function:: void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Computes the norm of ``op``.

    For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the norm
    of `a` as the determinant of this map.  Equivalently, if `\Sigma` generates
    `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of `a` is equal to
    `\prod_{i=0}^{d-1} \Sigma^i (a)`, where
    `d = \text{dim}_{\mathbf{F}_p}(\mathbf{F}_q)`.

    Algorithm selection is automatic depending on the input.

.. function:: void _fq_nmod_frobenius(mp_ptr *rop, mp_srcptr *op, slong len, slong e, const fq_nmod_ctx_t ctx)

    Sets ``(rop, 2d-1)`` to the image of ``(op, len)`` under the
    Frobenius operator raised to the e-th power, assuming that neither
    ``op`` nor ``e`` are zero.

.. function:: void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx)

    Evaluates the homomorphism `\Sigma^e` at ``op``.

    Recall that `\mathbf{F}_q / \mathbf{F}_p` is Galois with Galois group
    `\langle \sigma \rangle`, which is also isomorphic to
    `\mathbf{Z}/d\mathbf{Z}`, where
    `\sigma \in \operatorname{Gal}(\mathbf{F}_q/\mathbf{F}_p)` is the Frobenius element
    `\sigma \colon x \mapsto x^p`.

.. function:: int fq_nmod_multiplicative_order(fmpz_t ord, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Computes the order of ``op`` as an element of the
    multiplicative group of ``ctx``.
    
    Returns 0 if ``op`` is 0, otherwise it returns 1 if ``op``
    is a generator of the multiplicative group, and -1 if it is not.

    This function can also be used to check primitivity of a generator of
    a finite field whose defining polynomial is not primitive.

.. function:: int fq_nmod_is_primitive(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    Returns whether ``op`` is primitive, i.e., whether it is a
    generator of the multiplicative group of ``ctx``.


Bit packing
--------------------------------------------------------------------------------


.. function:: void fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, flint_bitcnt_t bit_size, const fq_nmod_ctx_t ctx)

    Packs ``op`` into bitfields of size ``bit_size``, writing the
    result to ``f``.

.. function:: void fq_nmod_bit_unpack(fq_nmod_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_nmod_ctx_t ctx)

    Unpacks into ``rop`` the element with coefficients packed into
    fields of size ``bit_size`` as represented by the integer
    ``f``.
