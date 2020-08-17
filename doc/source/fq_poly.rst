.. _fq-poly:

**fq_poly.h** -- univariate polynomials over finite fields
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fq_poly_struct

.. type:: fq_poly_t

    Description.


Memory management
--------------------------------------------------------------------------------


.. function:: void fq_poly_init(fq_poly_t poly, const fq_ctx_t ctx)

    Initialises ``poly`` for use, with context ctx, and setting its
    length to zero. A corresponding call to :func:`fq_poly_clear`
    must be made after finishing with the ``fq_poly_t`` to free the
    memory used by the polynomial.

.. function:: void fq_poly_init2(fq_poly_t poly, slong alloc, const fq_ctx_t ctx)

    Initialises ``poly`` with space for at least ``alloc``
    coefficients and sets the length to zero.  The allocated
    coefficients are all set to zero.  A corresponding call to
    :func:`fq_poly_clear` must be made after finishing with the
    ``fq_poly_t`` to free the memory used by the polynomial.

.. function:: void fq_poly_realloc(fq_poly_t poly, slong alloc, const fq_ctx_t ctx)

    Reallocates the given polynomial to have space for ``alloc``
    coefficients.  If ``alloc`` is zero the polynomial is cleared
    and then reinitialised.  If the current length is greater than
    ``alloc`` the polynomial is first truncated to length
    ``alloc``.

.. function:: void fq_poly_fit_length(fq_poly_t poly, slong len, const fq_ctx_t ctx)

    If ``len`` is greater than the number of coefficients currently
    allocated, then the polynomial is reallocated to have space for at
    least ``len`` coefficients.  No data is lost when calling this
    function.

    The function efficiently deals with the case where
    ``fit_length`` is called many times in small increments by at
    least doubling the number of allocated coefficients when length is
    larger than the number of coefficients currently allocated.

.. function:: void _fq_poly_set_length(fq_poly_t poly, slong newlen, const fq_ctx_t ctx)

    Sets the coefficients of ``poly`` beyond ``len`` to zero and
    sets the length of ``poly`` to ``len``.

.. function:: void fq_poly_clear(fq_poly_t poly, const fq_ctx_t ctx)

    Clears the given polynomial, releasing any memory used.  It must
    be reinitialised in order to be used again.

.. function:: void _fq_poly_normalise(fq_poly_t poly, const fq_ctx_t ctx)

    Sets the length of ``poly`` so that the top coefficient is
    non-zero.  If all coefficients are zero, the length is set to
    zero.  This function is mainly used internally, as all functions
    guarantee normalisation.

.. function:: void _fq_poly_normalise2(fq_struct *poly, slong *length, const fq_ctx_t ctx)

    Sets the length ``length`` of ``(poly,length)`` so that the
    top coefficient is non-zero. If all coefficients are zero, the
    length is set to zero. This function is mainly used internally, as
    all functions guarantee normalisation.

.. function:: void fq_poly_truncate(fq_poly_t poly, slong newlen, const fq_ctx_t ctx)

    Truncates the polynomial to length at most `n`.

.. function:: void fq_poly_set_trunc(fq_poly_t poly1, fq_poly_t poly2, slong newlen, const fq_ctx_t ctx)

    Sets ``poly1`` to ``poly2`` truncated to length `n`.

.. function:: void _fq_poly_reverse(fq_struct* output, const fq_struct* input, slong len, slong m, const fq_ctx_t ctx)

    Sets ``output`` to the reverse of ``input``, which is of
    length ``len``, but thinking of it as a polynomial of
    length ``m``, notionally zero-padded if necessary. The
    length ``m`` must be non-negative, but there are no other
    restrictions. The polynomial ``output`` must have space for
    ``m`` coefficients.

.. function:: void fq_poly_reverse(fq_poly_t output, const fq_poly_t input, slong m, const fq_ctx_t ctx)

    Sets ``output`` to the reverse of ``input``, thinking of it
    as a polynomial of length ``m``, notionally zero-padded if
    necessary).  The length ``m`` must be non-negative, but there
    are no other restrictions. The output polynomial will be set to
    length ``m`` and then normalised.


Polynomial parameters
--------------------------------------------------------------------------------


.. function:: long fq_poly_degree(fq_poly_t poly, const fq_ctx_t ctx)

    Returns the degree of the polynomial ``poly``.

.. function:: long fq_poly_length(fq_poly_t poly, const fq_ctx_t ctx)

    Returns the length of the polynomial ``poly``.

.. function:: fq_struct * fq_poly_lead(const fq_poly_t poly, const fq_ctx_t ctx)

    Returns a pointer to the leading coefficient of ``poly``, or
    ``NULL`` if ``poly`` is the zero polynomial.


Randomisation
--------------------------------------------------------------------------------


.. function:: void fq_poly_randtest(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

    Sets `f` to a random polynomial of length at most ``len``
    with entries in the field described by ``ctx``.

.. function:: void fq_poly_randtest_not_zero(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

    Same as ``fq_poly_randtest`` but guarantees that the polynomial
    is not zero.

.. function:: void fq_poly_randtest_monic(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

    Sets `f` to a random monic polynomial of length ``len`` with
    entries in the field described by ``ctx``.

.. function:: void fq_poly_randtest_irreducible(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

    Sets `f` to a random monic, irreducible polynomial of length
    ``len`` with entries in the field described by ``ctx``.


Assignment and basic manipulation
--------------------------------------------------------------------------------


.. function:: void _fq_poly_set(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, len``) to ``(op, len)``.

.. function:: void fq_poly_set(fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx)

    Sets the polynomial ``poly1`` to the polynomial ``poly2``.

.. function:: void fq_poly_set_fq(fq_poly_t poly, const fq_t c, const fq_ctx_t ctx)

    Sets the polynomial ``poly`` to ``c``.

.. function:: void fq_poly_set_fmpz_mod_poly(fq_poly_t rop, const fmpz_mod_poly_t op, fq_ctx_t ctx)

    Sets the polynomial ``rop`` to the polynomial ``op``

.. function:: void fq_poly_set_nmod_poly(fq_poly_t rop, const nmod_poly_t op, fq_ctx_t ctx)

    Sets the polynomial ``rop`` to the polynomial ``op``

.. function:: void fq_poly_swap(fq_poly_t op1, fq_poly_t op2, const fq_ctx_t ctx)

    Swaps the two polynomials ``op1`` and ``op2``.

.. function:: void _fq_poly_zero(fq_struct *rop, slong len, const fq_ctx_t ctx)

    Sets ``(rop, len)`` to the zero polynomial.

.. function:: void fq_poly_zero(fq_poly_t poly, const fq_ctx_t ctx)

    Sets ``poly`` to the zero polynomial.

.. function:: void fq_poly_one(fq_poly_t poly, const fq_ctx_t ctx)

    Sets ``poly`` to the constant polynomial `1`.

.. function:: void fq_poly_gen(fq_poly_t poly, const fq_ctx_t ctx)

    Sets ``poly`` to the polynomial `x`.

.. function:: void fq_poly_make_monic(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

     Sets ``rop`` to ``op``, normed to have leading coefficient 1.

.. function:: void _fq_poly_make_monic(fq_struct *rop, const fq_struct *op, slong length, const fq_ctx_t ctx)

     Sets ``rop`` to ``(op,length)``, normed to have leading coefficient 1.
     Assumes that ``rop`` has enough space for the polynomial, assumes that
     ``op`` is not zero (and thus has an invertible leading coefficient).


Getting and setting coefficients
--------------------------------------------------------------------------------


.. function:: void fq_poly_get_coeff(fq_t x, const fq_poly_t poly, slong n, const fq_ctx_t ctx)

    Sets `x` to the coefficient of `X^n` in ``poly``.

.. function:: void fq_poly_set_coeff(fq_poly_t poly, slong n, const fq_t x, const fq_ctx_t ctx)

    Sets the coefficient of `X^n` in ``poly`` to `x`.

.. function:: void fq_poly_set_coeff_fmpz(fq_poly_t poly, slong n, const fmpz_t x, const fq_ctx_t ctx)

    Sets the coefficient of `X^n` in the polynomial to `x`,
    assuming `n \geq 0`.


Comparison
--------------------------------------------------------------------------------


.. function:: int fq_poly_equal(const fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx)

    Returns nonzero if the two polynomials ``poly1`` and ``poly2``
    are equal, otherwise returns zero.

.. function:: int fq_poly_equal_trunc(const fq_poly_t poly1, const fq_poly_t poly2, slong n, const fq_ctx_t ctx)

    Notionally truncate ``poly1`` and ``poly2`` to length `n` and
    return nonzero if they are equal, otherwise return zero.

.. function:: int fq_poly_is_zero(const fq_poly_t poly, const fq_ctx_t ctx)

    Returns whether the polynomial ``poly`` is the zero polynomial.

.. function:: int fq_poly_is_one(const fq_poly_t op)

    Returns whether the polynomial ``poly`` is equal
    to the constant polynomial `1`.

.. function:: int fq_poly_is_gen(const fq_poly_t op, const fq_ctx_t ctx)

    Returns whether the polynomial ``poly`` is equal
    to the polynomial `x`.

.. function:: int fq_poly_is_unit(const fq_poly_t op, const fq_ctx_t ctx)

    Returns whether the polynomial ``poly`` is a unit in the polynomial
    ring `\mathbf{F}_q[X]`, i.e. if it has degree `0` and is non-zero.

.. function:: int fq_poly_equal_fq(const fq_poly_t poly, const fq_t c, const fq_ctx_t ctx)

    Returns whether the polynomial ``poly`` is equal the (constant)
    `\mathbf{F}_q` element ``c``


Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void _fq_poly_add(fq_struct *res, const fq_struct *poly1, slong len1, const fq_struct *poly2, slong len2, const fq_ctx_t ctx)

    Sets ``res`` to the sum of ``(poly1,len1)`` and ``(poly2,len2)``.

.. function:: void fq_poly_add(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx)

    Sets ``res`` to the sum of ``poly1`` and ``poly2``.

.. function:: void fq_poly_add_si(fq_poly_t res, const fq_poly_t poly1, slong c, const fq_ctx_t ctx)

    Sets ``res`` to the sum of ``poly1`` and ``c``.

.. function:: void fq_poly_add_series(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, slong n, const fq_ctx_t ctx)

    Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    ``res`` to the sum.

.. function:: void _fq_poly_sub(fq_struct *res, const fq_struct *poly1, slong len1, const fq_struct *poly2, slong len2, const fq_ctx_t ctx)

    Sets ``res`` to the difference of ``(poly1,len1)`` and ``(poly2,len2)``.

.. function:: void fq_poly_sub(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx)

    Sets ``res`` to the difference of ``poly1`` and ``poly2``.

.. function:: void fq_poly_sub_series(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, slong n, const fq_ctx_t ctx)

    Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    ``res`` to the difference.

.. function:: void _fq_poly_neg(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``rop`` to the additive inverse of ``(poly,len)``.

.. function:: void fq_poly_neg(fq_poly_t res, const fq_poly_t poly, const fq_ctx_t ctx)

    Sets ``res`` to the additive inverse of ``poly``.


Scalar multiplication and division
--------------------------------------------------------------------------------


.. function:: void _fq_poly_scalar_mul_fq(fq_struct *rop, const fq_struct *op, slong len, const fq_t x, const fq_ctx_t ctx)

    Sets ``(rop,len)`` to the product of ``(op,len)`` by the
    scalar ``x``, in the context defined by ``ctx``.

.. function:: void fq_poly_scalar_mul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op`` by the scalar ``x``, in the context
    defined by ``ctx``.

.. function:: void _fq_poly_scalar_addmul_fq(fq_struct *rop, const fq_struct *op, slong len, const fq_t x, const fq_ctx_t ctx)

    Adds to ``(rop,len)`` the product of ``(op,len)`` by the
    scalar ``x``, in the context defined by ``ctx``.
    In particular, assumes the same length for ``op`` and
    ``rop``.

.. function:: void fq_poly_scalar_addmul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

    Adds to ``rop`` the product of ``op`` by the
    scalar ``x``, in the context defined by ``ctx``.

.. function:: void _fq_poly_scalar_submul_fq(fq_struct *rop, const fq_struct *op, slong len, const fq_t x, const fq_ctx_t ctx)

    Subtracts from ``(rop,len)`` the product of ``(op,len)`` by the
    scalar ``x``, in the context defined by ``ctx``.
    In particular, assumes the same length for ``op`` and
    ``rop``.

.. function:: void fq_poly_scalar_submul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

    Subtracts from ``rop`` the product of ``op`` by the
    scalar ``x``, in the context defined by ``ctx``.

.. function:: void _fq_poly_scalar_div_fq(fq_struct *rop, const fq_struct *op, slong len, const fq_t x, const fq_ctx_t ctx)

    Sets ``(rop,len)`` to the quotient of ``(op,len)`` by the
    scalar ``x``, in the context defined by ``ctx``. An exception is raised
    if ``x`` is zero.

.. function:: void fq_poly_scalar_div_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)                                                 

    Sets ``rop`` to the quotient of ``op`` by the scalar ``x``, in the context
    defined by ``ctx``. An exception is raised if ``x`` is zero.

Multiplication
--------------------------------------------------------------------------------


.. function:: void _fq_poly_mul_classical(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    and ``(op2, len2)``, assuming that ``len1`` is at least ``len2``
    and neither is zero.

    Permits zero padding.  Does not support aliasing of ``rop``
    with either ``op1`` or ``op2``.

.. function:: void fq_poly_mul_classical(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``
    using classical polynomial multiplication.

.. function:: void _fq_poly_mul_reorder(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    and ``(op2, len2)``, assuming that ``len1`` and ``len2`` are
    non-zero.

    Permits zero padding.  Supports aliasing.

.. function:: void fq_poly_mul_reorder(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``,
    reordering the two indeterminates `X` and `Y` when viewing
    the polynomials as elements of `\mathbf{F}_p[X,Y]`.

    Suppose `\mathbf{F}_q = \mathbf{F}_p[X]/ (f(X))` and recall
    that elements of `\mathbf{F}_q` are internally represented
    by elements of type ``fmpz_poly``.  For small degree extensions
    but polynomials in `\mathbf{F}_q[Y]` of large degree `n`, we
    change the representation to

    .. math ::


        \begin{split}
        g(Y) & = \sum_{i=0}^{n} a_i(X) Y^i \\
             & = \sum_{j=0}^{d} \sum_{i=0}^{n} \text{Coeff}(a_i(X), j) Y^i.
        \end{split}


    This allows us to use a poor algorithm (such as classical multiplication)
    in the `X`-direction and leverage the existing fast integer
    multiplication routines in the `Y`-direction where the polynomial
    degree `n` is large.

.. function:: void _fq_poly_mul_univariate(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    and ``(op2, len2)``.

    Permits zero padding and places no assumptions on the
    lengths ``len1`` and ``len2``.  Supports aliasing.

.. function:: void fq_poly_mul_univariate(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``
    using a bivariate to univariate transformation and reducing
    this problem to multiplying two univariate polynomials.

.. function:: void _fq_poly_mul_KS(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    and ``(op2, len2)``.

    Permits zero padding and places no assumptions on the
    lengths ``len1`` and ``len2``.  Supports aliasing.

.. function:: void fq_poly_mul_KS(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``
    using Kronecker substitution, that is, by encoding each
    coefficient in `\mathbf{F}_{q}` as an integer and reducing
    this problem to multiplying two polynomials over the integers.

.. function:: void _fq_poly_mul(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``(rop, len1 + len2 - 1)`` to the product of ``(op1, len1)``
    and ``(op2, len2)``, choosing an appropriate algorithm.

    Permits zero padding.  Does not support aliasing.

.. function:: void fq_poly_mul(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``op1`` and ``op2``,
    choosing an appropriate algorithm.

.. function:: void _fq_poly_mullow_classical(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, slong n, const fq_ctx_t ctx)

    Sets ``(rop, n)`` to the first `n` coefficients of ``(op1, len1)``
    multiplied by ``(op2, len2)``.

    Assumes ``0 < n <= len1 + len2 - 1``.  Assumes neither ``len1`` nor
    ``len2`` is zero.

.. function:: void fq_poly_mullow_classical(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to the product of ``poly1`` and ``poly2``, computed
    using the classical or schoolbook method.

.. function:: void _fq_poly_mullow_univariate(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, slong n, const fq_ctx_t ctx)

    Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    ``(op1, len1)`` and ``(op2, len2)``, computed using a
    bivariate to univariate transformation.

    Assumes that ``len1`` and ``len2`` are positive, but does allow
    for the polynomials to be zero-padded.  The polynomials may be zero,
    too.  Assumes `n` is positive.  Supports aliasing between ``res``,
    ``poly1`` and ``poly2``.

.. function:: void fq_poly_mullow_univariate(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to the lowest `n` coefficients of the product of
    ``op1`` and ``op2``, computed using a bivariate to univariate
    transformation.

.. function:: void _fq_poly_mullow_KS(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, slong n, const fq_ctx_t ctx)

    Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    ``(op1, len1)`` and ``(op2, len2)``.

    Assumes that ``len1`` and ``len2`` are positive, but does allow
    for the polynomials to be zero-padded.  The polynomials may be zero,
    too.  Assumes `n` is positive.  Supports aliasing between ``rop``,
    ``op1`` and ``op2``.

.. function:: void fq_poly_mullow_KS(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to the lowest `n` coefficients of the product of
    ``op1`` and ``op2``.

.. function:: void _fq_poly_mullow(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, slong n, const fq_ctx_t ctx)

    Sets ``(rop, n)`` to the lowest `n` coefficients of the product of
    ``(op1, len1)`` and ``(op2, len2)``.

    Assumes ``0 < n <= len1 + len2 - 1``.  Allows for zero-padding in
    the inputs.  Does not support aliasing between the inputs and the output.

.. function:: void fq_poly_mullow(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to the lowest `n` coefficients of the product of
    ``op1`` and ``op2``.

.. function:: void _fq_poly_mulhigh_classical(fq_struct *res, const fq_struct *poly1, slong len1, const fq_struct *poly2, slong len2, slong start, const fq_ctx_t ctx)

    Computes the product of ``(poly1, len1)`` and ``(poly2, len2)``
    and writes the coefficients from ``start`` onwards into the high
    coefficients of ``res``, the remaining coefficients being arbitrary
    but reduced.  Assumes that ``len1 >= len2 > 0``. Aliasing of inputs
    and output is not permitted.  Algorithm is classical multiplication.

.. function:: void fq_poly_mulhigh_classical(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, slong start, const fq_ctx_t ctx)

    Computes the product of ``poly1`` and ``poly2`` and writes the
    coefficients from ``start`` onwards into the high coefficients of
    ``res``, the remaining coefficients being arbitrary but reduced.
    Algorithm is classical multiplication.

.. function:: void _fq_poly_mulhigh(fq_struct *res, const fq_struct *poly1, slong len1, const fq_struct *poly2, slong len2, slong start, const fq_ctx_t ctx)

    Computes the product of ``(poly1, len1)`` and ``(poly2, len2)``
    and writes the coefficients from ``start`` onwards into the high
    coefficients of ``res``, the remaining coefficients being arbitrary
    but reduced.  Assumes that ``len1 >= len2 > 0``. Aliasing of inputs
    and output is not permitted.

.. function:: void fq_poly_mulhigh(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, slong start, const fq_ctx_t ctx)

    Computes the product of ``poly1`` and ``poly2`` and writes the
    coefficients from ``start`` onwards into the high coefficients of
    ``res``, the remaining coefficients being arbitrary but reduced.

.. function:: void _fq_poly_mulmod(fq_struct* res, const fq_struct* poly1, slong len1, const fq_struct* poly2, slong len2, const fq_struct* f, slong lenf, const fq_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

    It is required that ``len1 + len2 - lenf > 0``, which is
    equivalent to requiring that the result will actually be
    reduced. Otherwise, simply use ``_fq_poly_mul`` instead.

    Aliasing of ``f`` and ``res`` is not permitted.

.. function:: void fq_poly_mulmod(fq_poly_t res,const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t f, const fq_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

.. function:: void _fq_poly_mulmod_preinv(fq_struct* res, const fq_struct* poly1, slong len1, const fq_struct* poly2, slong len2, const fq_struct* f, slong lenf, const fq_struct* finv, slong lenfinv, const fq_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

    It is required that ``finv`` is the inverse of the reverse of
    ``f`` mod ``x^lenf``.

    Aliasing of ``res`` with any of the inputs is not permitted.

.. function:: void fq_poly_mulmod_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``. ``finv``
    is the inverse of the reverse of ``f``.


Squaring
--------------------------------------------------------------------------------


.. function:: void _fq_poly_sqr_classical(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, 2*len - 1)`` to the square of ``(op, len)``,
    assuming that ``(op,len)`` is not zero and using classical
    polynomial multiplication.

    Permits zero padding.  Does not support aliasing of ``rop``
    with either ``op1`` or ``op2``.

.. function:: void fq_poly_sqr_classical(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

    Sets ``rop`` to the square of ``op`` using classical
     polynomial multiplication.


.. function:: void _fq_poly_sqr_reorder(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, 2*len- 1)`` to the square of ``(op, len)``,
    assuming that ``len`` is not zero reordering the two indeterminates
    `X` and `Y` when viewing the polynomials as elements of `\mathbf{F}_p[X,Y]`.

    Permits zero padding.  Supports aliasing.

.. function:: void fq_poly_sqr_reorder(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

    Sets ``rop`` to the square of ``op``,
    assuming that ``len`` is not zero reordering the two indeterminates
    `X` and `Y` when viewing the polynomials as elements of `\mathbf{F}_p[X,Y]`.
    See ``fq_poly_mul_reorder``.


.. function:: void _fq_poly_sqr_KS(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, 2*len - 1)`` to the square of ``(op, len)``.

    Permits zero padding and places no assumptions on the
    lengths ``len1`` and ``len2``.  Supports aliasing.

.. function:: void fq_poly_sqr_KS(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

    Sets ``rop`` to the square ``op`` using Kronecker substitution,
    that is, by encoding each coefficient in `\mathbf{F}_{q}` as an integer
    and reducing this problem to multiplying two polynomials over the integers.

.. function:: void _fq_poly_sqr(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, 2* len - 1)`` to the square of ``(op, len)``,
    choosing an appropriate algorithm.

    Permits zero padding.  Does not support aliasing.

.. function:: void fq_poly_sqr(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

    Sets ``rop`` to the square of ``op``,
    choosing an appropriate algorithm.



Powering
--------------------------------------------------------------------------------


.. function:: void _fq_poly_pow(fq_struct *rop, const fq_struct *op, slong len, ulong e, const fq_ctx_t ctx)

    Sets ``rop = op^e``, assuming that ``e, len > 0`` and that
    ``rop`` has space for ``e*(len - 1) + 1`` coefficients.  Does
    not support aliasing.

.. function:: void fq_poly_pow(fq_poly_t rop, const fq_poly_t op, ulong e, const fq_ctx_t ctx)

    Computes ``rop = op^e``.  If `e` is zero, returns one,
    so that in particular ``0^0 = 1``.

.. function:: void _fq_poly_powmod_ui_binexp(fq_struct* res, const fq_struct* poly, ulong e, const fq_struct* f, slong lenf, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e > 0``.

    We require ``lenf > 1``. It is assumed that ``poly`` is
    already reduced modulo ``f`` and zero-padded as necessary to
    have length exactly ``lenf - 1``. The output ``res`` must
    have room for ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_ui_binexp(fq_poly_t res, const fq_poly_t poly, ulong e, const fq_poly_t f, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e >= 0``.

.. function:: void _fq_poly_powmod_ui_binexp_preinv(fq_struct* res, const fq_struct* poly, ulong e, const fq_struct* f, slong lenf, const fq_struct* finv, slong lenfinv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e > 0``.
    We require ``finv`` to be the inverse of the reverse of
    ``f``.

    We require ``lenf > 1``. It is assumed that ``poly`` is
    already reduced modulo ``f`` and zero-padded as necessary to
    have length exactly ``lenf - 1``. The output ``res`` must
    have room for ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_ui_binexp_preinv(fq_poly_t res, const fq_poly_t poly, ulong e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e >= 0``.
    We require ``finv`` to be the inverse of the reverse of
    ``f``.

.. function:: void _fq_poly_powmod_fmpz_binexp(fq_struct* res, const fq_struct* poly, fmpz_t e, const fq_struct* f, slong lenf, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e > 0``.

    We require ``lenf > 1``. It is assumed that ``poly`` is
    already reduced modulo ``f`` and zero-padded as necessary to
    have length exactly ``lenf - 1``. The output ``res`` must
    have room for ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_fmpz_binexp(fq_poly_t res, const fq_poly_t poly, fmpz_t e, const fq_poly_t f, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e >= 0``.

.. function:: void _fq_poly_powmod_fmpz_binexp_preinv(fq_struct* res, const fq_struct* poly, fmpz_t e, const fq_struct* f, slong lenf, const fq_struct* finv, slong lenfinv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e > 0``.
    We require ``finv`` to be the inverse of the reverse of
    ``f``.

    We require ``lenf > 1``. It is assumed that ``poly`` is
    already reduced modulo ``f`` and zero-padded as necessary to
    have length exactly ``lenf - 1``. The output ``res`` must
    have room for ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_fmpz_binexp_preinv(fq_poly_t res, const fq_poly_t poly, fmpz_t e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using binary exponentiation. We require ``e >= 0``.
    We require ``finv`` to be the inverse of the reverse of
    ``f``.

.. function:: void _fq_poly_powmod_fmpz_sliding_preinv(fq_struct* res, const fq_struct* poly, fmpz_t e, ulong k, const fq_struct* f, slong lenf, const fq_struct* finv, slong lenfinv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using sliding-window exponentiation with window size
    ``k``. We require ``e > 0``.  We require ``finv`` to be
    the inverse of the reverse of ``f``. If ``k`` is set to
    zero, then an "optimum" size will be selected automatically base
    on ``e``.

    We require ``lenf > 1``. It is assumed that ``poly`` is
    already reduced modulo ``f`` and zero-padded as necessary to
    have length exactly ``lenf - 1``. The output ``res`` must
    have room for ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_fmpz_sliding_preinv(fq_poly_t res, const fq_poly_t poly, fmpz_t e, ulong k, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    ``f``, using sliding-window exponentiation with window size
    ``k``. We require ``e >= 0``.  We require ``finv`` to be
    the inverse of the reverse of ``f``.  If ``k`` is set to
    zero, then an "optimum" size will be selected automatically base
    on ``e``.

.. function:: void _fq_poly_powmod_x_fmpz_preinv(fq_struct * res, const fmpz_t e, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)

    Sets ``res`` to ``x`` raised to the power ``e`` modulo ``f``,
    using sliding window exponentiation. We require ``e > 0``.
    We require ``finv`` to be the inverse of the reverse of ``f``.

    We require ``lenf > 2``. The output ``res`` must have room for
    ``lenf - 1`` coefficients.

.. function:: void fq_poly_powmod_x_fmpz_preinv(fq_poly_t res, const fmpz_t e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

    Sets ``res`` to ``x`` raised to the power ``e``
    modulo ``f``, using sliding window exponentiation. We require
    ``e >= 0``. We require ``finv`` to be the inverse of the reverse of
    ``f``.

.. function:: void _fq_poly_pow_trunc_binexp(fq_struct * res, const fq_struct * poly, ulong e, slong trunc, const fq_ctx_t ctx)

    Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    (assumed to be zero padded if necessary to length ``trunc``) to
    the power ``e``. This is equivalent to doing a powering followed
    by a truncation. We require that ``res`` has enough space for
    ``trunc`` coefficients, that ``trunc > 0`` and that
    ``e > 1``. Aliasing is not permitted. Uses the binary
    exponentiation method.

.. function:: void fq_poly_pow_trunc_binexp(fq_poly_t res, const fq_poly_t poly, ulong e, slong trunc, const fq_ctx_t ctx)

    Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    to the power ``e``. This is equivalent to doing a powering
    followed by a truncation. Uses the binary exponentiation method.

.. function:: void _fq_poly_pow_trunc(fq_struct * res, const fq_struct * poly, ulong e, slong trunc, const fq_ctx_t mod)

    Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    (assumed to be zero padded if necessary to length ``trunc``) to
    the power ``e``. This is equivalent to doing a powering followed
    by a truncation. We require that ``res`` has enough space for
    ``trunc`` coefficients, that ``trunc > 0`` and that
    ``e > 1``. Aliasing is not permitted.

.. function:: void fq_poly_pow_trunc(fq_poly_t res, const fq_poly_t poly, ulong e, slong trunc, fq_ctx_t ctx)

    Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    to the power ``e``. This is equivalent to doing a powering
    followed by a truncation.


Shifting
--------------------------------------------------------------------------------


.. function:: void _fq_poly_shift_left(fq_struct *rop, const fq_struct *op, slong len, slong n, const fq_ctx_t ctx)

    Sets ``(rop, len + n)`` to ``(op, len)`` shifted left by
    `n` coefficients.

    Inserts zero coefficients at the lower end.  Assumes that
    ``len`` and `n` are positive, and that ``rop`` fits
    ``len + n`` elements.  Supports aliasing between ``rop`` and
    ``op``.

.. function:: void fq_poly_shift_left(fq_poly_t rop, const fq_poly_t op, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to ``op`` shifted left by `n` coeffs.  Zero
    coefficients are inserted.

.. function:: void _fq_poly_shift_right(fq_struct *rop, const fq_struct *op, slong len, slong n, const fq_ctx_t ctx)

    Sets ``(rop, len - n)`` to ``(op, len)`` shifted right by
    `n` coefficients.

    Assumes that ``len`` and `n` are positive, that ``len > n``,
    and that ``rop`` fits ``len - n`` elements.  Supports
    aliasing between ``rop`` and ``op``, although in this case
    the top coefficients of ``op`` are not set to zero.

.. function:: void fq_poly_shift_right(fq_poly_t rop, const fq_poly_t op, slong n, const fq_ctx_t ctx)

    Sets ``rop`` to ``op`` shifted right by `n` coefficients.
    If `n` is equal to or greater than the current length of
    ``op``, ``rop`` is set to the zero polynomial.


Norms
--------------------------------------------------------------------------------


.. function:: long _fq_poly_hamming_weight(const fq_poly *op, slong len, const fq_ctx_t ctx)

    Returns the number of non-zero entries in ``(op, len)``.

.. function:: long fq_poly_hamming_weight(const fq_poly_t op, const fq_ctx_t ctx)

    Returns the number of non-zero entries in the polynomial ``op``.


Euclidean division
--------------------------------------------------------------------------------


.. function:: void _fq_poly_divrem_basecase(fq_struct *Q, fq_struct *R, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Computes ``(Q, lenA - lenB + 1)``, ``(R, lenA)`` such that
    `A = B Q + R` with `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible
    and that ``invB`` is its inverse.

    Assumes that `\operatorname{len}(A), \operatorname{len}(B) > 0`.  Allows zero-padding in
    ``(A, lenA)``.  `R` and `A` may be aliased, but apart from
    this no aliasing of input and output operands is allowed.

.. function:: void fq_poly_divrem_basecase(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Computes `Q`, `R` such that `A = B Q + R` with
    `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible.  This can
    be taken for granted the context is for a finite field, that is, when
    `p` is prime and `f(X)` is irreducible.

.. function:: void _fq_poly_divrem(fq_struct *Q, fq_struct *R, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Computes ``(Q, lenA - lenB + 1)``, ``(R, lenA)`` such that
    `A = B Q + R` with `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible
    and that ``invB`` is its inverse.

    Assumes that `\operatorname{len}(A), \operatorname{len}(B) > 0`.  Allows zero-padding in
    ``(A, lenA)``.  `R` and `A` may be aliased, but apart from
    this no aliasing of input and output operands is allowed.

.. function:: void fq_poly_divrem(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Computes `Q`, `R` such that `A = B Q + R` with
    `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible.  This can
    be taken for granted the context is for a finite field, that is, when
    `p` is prime and `f(X)` is irreducible.

.. function:: void fq_poly_divrem_f(fq_t f, fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Either finds a non-trivial factor `f` of the modulus of
    ``ctx``, or computes `Q`, `R` such that `A = B Q + R` and
    `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    If the leading coefficient of `B` is invertible, the division with
    remainder operation is carried out, `Q` and `R` are computed
    correctly, and `f` is set to `1`.  Otherwise, `f` is set to a
    non-trivial factor of the modulus and `Q` and `R` are not touched.

    Assumes that `B` is non-zero.

.. function:: void _fq_poly_rem(fq_struct *R, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Sets ``R`` to the remainder of the division of ``(A,lenA)`` by
    ``(B,lenB)``. Assumes that the leading coefficient of ``(B,lenB)``
    is invertible and that ``invB`` is its inverse.

.. function:: void fq_poly_rem(fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Sets ``R`` to the remainder of the division of ``A`` by
    ``B`` in the context described by ``ctx``.

.. function:: void _fq_poly_div_basecase(fq_struct *Q, fq_struct *R, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Notationally, computes `Q`, `R` such that `A = B Q + R` with `0
    \leq \operatorname{len}(R) < \operatorname{len}(B)` but only sets ``(Q, lenA - lenB + 1)``.

    Requires temporary space ``(R, lenA)``.  If ``R`` is
    ``NULL``, then the temporary space will be allocated.  Allows
    aliasing only between `A` and `R`.  Allows zero-padding in `A` but
    not in `B`.  Assumes that the leading coefficient of `B` is a
    unit.

.. function:: void fq_poly_div_basecase(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Notionally finds polynomials `Q` and `R` such that `A = B Q + R` with
    `\operatorname{len}(R) < \operatorname{len}(B)`, but returns only ``Q``. If `\operatorname{len}(B) = 0` an
    exception is raised.

.. function:: void _fq_poly_divrem_divconquer_recursive(fq_struct * Q, fq_struct * BQ, fq_struct * W, const fq_struct * A, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Computes ``(Q, lenB)``, ``(BQ, 2 lenB - 1)`` such that
    `BQ = B \times Q` and `A = B Q + R` where `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible and that
    ``invB`` is the inverse.

    Assumes `\operatorname{len}(B) > 0`.  Allows zero-padding in ``(A, lenA)``.  Requires
    a temporary array ``(W, 2 lenB - 1)``.  No aliasing of input and output
    operands is allowed.

    This function does not read the bottom `\operatorname{len}(B) - 1` coefficients from
    `A`, which means that they might not even need to exist in allocated
    memory.

.. function:: void _fq_poly_divrem_divconquer(fq_struct * Q, fq_struct * R, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Computes ``(Q, lenA - lenB + 1)``, ``(R, lenA)`` such that
    `A = B Q + R` and `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that the leading coefficient of `B` is invertible and that
    ``invB`` is the inverse.

    Assumes `\operatorname{len}(A) \geq \operatorname{len}(B) > 0`.  Allows zero-padding in
    ``(A, lenA)``.  No aliasing of input and output operands is
    allowed.

.. function:: void fq_poly_divrem_divconquer(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Computes `Q`, `R` such that `A = B Q + R` and `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.

    Assumes that `B` is non-zero and that the leading coefficient of
    `B` is invertible.

.. function:: void _fq_poly_div_newton_n_preinv(fq_struct* Q, const fq_struct* A, slong lenA, const fq_struct* B, slong lenB, const fq_struct* Binv, slong lenBinv, const fq_struct ctx_t)

    Notionally computes polynomials `Q` and `R` such that `A = BQ + R` with
    `\operatorname{len}(R)` less than ``lenB``, where ``A`` is of length ``lenA``
    and ``B`` is of length ``lenB``, but return only `Q`.

    We require that `Q` have space for ``lenA - lenB + 1`` coefficients
    and assume that the leading coefficient of `B` is a unit. Furthermore, we
    assume that `Binv` is the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.

    The algorithm used is to reverse the polynomials and divide the
    resulting power series, then reverse the result.

.. function:: void fq_poly_div_newton_n_preinv(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_poly_t Binv, const fq_ctx_t ctx)

    Notionally computes `Q` and `R` such that `A = BQ + R` with
    `\operatorname{len}(R) < \operatorname{len}(B)`, but returns only `Q`.

    We assume that the leading coefficient of `B` is a unit and that `Binv` is
    the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.

    It is required that the length of `A` is less than or equal to
    2*the length of `B` - 2.

    The algorithm used is to reverse the polynomials and divide the
    resulting power series, then reverse the result.

.. function:: void _fq_poly_divrem_newton_n_preinv(fq_struct* Q, fq_struct* R, const fq_struct* A, slong lenA, const fq_struct* B, slong lenB, const fq_struct* Binv, slong lenBinv, const fq_ctx_t ctx)

    Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R)` less
    than ``lenB``, where `A` is of length ``lenA`` and `B` is of
    length ``lenB``. We require that `Q` have space for
    ``lenA - lenB + 1`` coefficients. Furthermore, we assume that `Binv` is
    the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`. The algorithm
    used is to call :func:`div_newton_n_preinv` and then multiply out
    and compute the remainder.

.. function:: void fq_poly_divrem_newton_preinv(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_poly_t Binv, const fq_ctx_t ctx)

    Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R) <
    \operatorname{len}(B)`.  We assume `Binv` is the inverse of the reverse of `B`
    mod `x^{\operatorname{len}(B)}`.

    It is required that the length of `A` is less than or equal to
    2*the length of `B` - 2.

    The algorithm used is to call :func:`div_newton_n` and then
    multiply out and compute the remainder.

.. function:: void _fq_poly_inv_series_newton(fq_struct* Qinv, const fq_struct* Q, slong n, const fq_ctx_t ctx)

    Given ``Q`` of length ``n`` whose constant coefficient is
    invertible modulo the given modulus, find a polynomial ``Qinv``
    of length ``n`` such that ``Q * Qinv`` is ``1`` modulo
    `x^n`. Requires ``n > 0``.  This function can be viewed as
    inverting a power series via Newton iteration.

.. function:: void fq_poly_inv_series_newton(fq_poly_t Qinv, const fq_poly_t Q, slong n, const fq_ctx_t ctx)

    Given ``Q`` find ``Qinv`` such that ``Q * Qinv`` is
    ``1`` modulo `x^n`. The constant coefficient of ``Q`` must
    be invertible modulo the modulus of ``Q``. An exception is
    raised if this is not the case or if ``n = 0``. This function
    can be viewed as inverting a power series via Newton iteration.

.. function:: void _fq_poly_inv_series(fq_struct* Qinv, const fq_struct* Q, slong n, const fq_ctx_t ctx)

    Given ``Q`` of length ``n`` whose constant coefficient is
    invertible modulo the given modulus, find a polynomial ``Qinv``
    of length ``n`` such that ``Q * Qinv`` is ``1`` modulo
    `x^n`. Requires ``n > 0``.

.. function:: void fq_poly_inv_series(fq_poly_t Qinv, const fq_poly_t Q, slong n, const fq_ctx_t ctx)

    Given ``Q`` find ``Qinv`` such that ``Q * Qinv`` is
    ``1`` modulo `x^n`. The constant coefficient of ``Q`` must
    be invertible modulo the modulus of ``Q``. An exception is
    raised if this is not the case or if ``n = 0``.

.. function:: void _fq_poly_div_series(fmpz * Q, const fmpz * A, slong Alen, const fmpz * B, slong Blen, slong n, fq_ctx_t ctx)

    Set ``(Q, n)`` to the quotient of the series ``(A, Alen``) and
    ``(B, Blen)`` assuming ``Alen, Blen <= n``. We assume the bottom
    coefficient of ``B`` is invertible.

.. function:: void fq_poly_div_series(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, slong n, fq_ctx_t ctx)

    Set `Q` to the quotient of the series `A` by `B`, thinking of the series as
    though they were of length `n`. We assume that the bottom coefficient of
    `B` is invertible.


Greatest common divisor
--------------------------------------------------------------------------------


.. function:: void fq_poly_gcd(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the greatest common divisor of ``op1`` and
    ``op2``, using the either the Euclidean or HGCD algorithm. The
    GCD of zero polynomials is defined to be zero, whereas the GCD of
    the zero polynomial and some other polynomial `P` is defined to be
    `P`. Except in the case where the GCD is zero, the GCD `G` is made
    monic.

.. function:: long _fq_poly_gcd(fq_struct* G,const fq_struct* A, slong lenA, const fq_struct* B, slong lenB, const fq_ctx_t ctx)

    Computes the GCD of `A` of length ``lenA`` and `B` of length
    ``lenB``, where ``lenA >= lenB > 0`` and sets `G` to it. The
    length of the GCD `G` is returned by the function. No attempt is
    made to make the GCD monic. It is required that `G` have space for
    ``lenB`` coefficients.

.. function:: void fq_poly_gcd_euclidean(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the greatest common divisor of ``op1`` and
    ``op2``, using the Euclidean algorithm. The GCD of zero
    polynomials is defined to be zero, whereas the GCD of the zero
    polynomial and some other polynomial `P` is defined to be
    `P`. Except in the case where the GCD is zero, the GCD `G` is made
    monic.

.. function:: long _fq_poly_gcd_euclidean(fq_struct* G, const fq_struct* A, slong lenA, const fq_struct* B, slong lenB, const fq_ctx_t ctx)

    Computes the GCD of `A` of length ``lenA`` and `B` of length
    ``lenB``, where ``lenA >= lenB > 0`` and sets `G` to it. The
    length of the GCD `G` is returned by the function. No attempt is
    made to make the GCD monic. It is required that `G` have space for
    ``lenB`` coefficients.

.. function:: slong _fq_poly_hgcd(fq_struct **M, slong *lenM, fq_struct *A, slong *lenA, fq_struct *B, slong *lenB, const fq_struct * a, slong lena, const fq_struct *b, slong lenb, const fq_ctx_t ctx)

    Computes the HGCD of `a` and `b`, that is, a matrix `M`, a sign `\sigma`
    and two polynomials `A` and `B` such that

    .. math ::


        (A,B)^t = \sigma M^{-1} (a,b)^t.



    Assumes that `\operatorname{len}(a) > \operatorname{len}(b) > 0`.

    Assumes that `A` and `B` have space of size at least `\operatorname{len}(a)`
    and `\operatorname{len}(b)`, respectively.  On exit, ``*lenA`` and ``*lenB``
    will contain the correct lengths of `A` and `B`.

    Assumes that ``M[0]``, ``M[1]``, ``M[2]``, and ``M[3]``
    each point to a vector of size at least `\operatorname{len}(a)`.

.. function:: void fq_poly_gcd_hgcd(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the greatest common divisor of ``op1`` and
    ``op2``, using the HGCD algorithm. The GCD of zero
    polynomials is defined to be zero, whereas the GCD of the zero
    polynomial and some other polynomial `P` is defined to be
    `P`. Except in the case where the GCD is zero, the GCD `G` is made
    monic.

.. function:: long _fq_poly_gcd_hgcd(fq_struct* G, const fq_struct* A, slong lenA, const fq_struct* B, slong lenB, const fq_ctx_t ctx)

    Computes the GCD of `A` of length ``lenA`` and `B` of length
    ``lenB`` using the HGCD algorithm, where
    ``lenA >= lenB > 0`` and sets `G` to it. The length of the GCD
    `G` is returned by the function. No attempt is made to make the
    GCD monic. It is required that `G` have space for ``lenB``
    coefficients.

.. function:: slong _fq_poly_gcd_euclidean_f(fq_t f, fq_struct *G, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_ctx_t ctx)

    Either sets `f = 1` and `G` to the greatest common divisor of
    `(A,\operatorname{len}(A))` and `(B, \operatorname{len}(B))` and returns its length, or sets
    `f` to a non-trivial factor of the modulus of ``ctx`` and leaves
    the contents of the vector `(G, lenB)` undefined.

    Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) > 0` and that the vector `G`
    has space for sufficiently many coefficients.

.. function:: void fq_poly_gcd_euclidean_f(fq_t f, fq_poly_t G, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Either sets `f = 1` and `G` to the greatest common divisor of `A`
    and `B` or sets `f` to a factor of the modulus of ``ctx``.

.. function:: slong _fq_poly_xgcd_euclidean(fq_struct *G, fq_struct *S, fq_struct *T, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fmpz_t invB, const fq_ctx_t ctx)

    Computes the GCD of `A` and `B` together with cofactors `S` and `T`
    such that `S A + T B = G`.  Returns the length of `G`.

    Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) \geq 1` and
    `(\operatorname{len}(A),\operatorname{len}(B)) \neq (1,1)`.

    No attempt is made to make the GCD monic.

    Requires that `G` have space for `\operatorname{len}(B)` coefficients.  Writes
    `\operatorname{len}(B)-1` and `\operatorname{len}(A)-1` coefficients to `S` and `T`, respectively.
    Note that, in fact, `\operatorname{len}(S) \leq \max(\operatorname{len}(B) - \operatorname{len}(G), 1)` and
    `\operatorname{len}(T) \leq \max(\operatorname{len}(A) - \operatorname{len}(G), 1)`.

    No aliasing of input and output operands is permitted.

.. function:: void fq_poly_xgcd_euclidean(fq_poly_t G, fq_poly_t S, fq_poly_t T, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Computes the GCD of `A` and `B`. The GCD of zero polynomials is
    defined to be zero, whereas the GCD of the zero polynomial and some other
    polynomial `P` is defined to be `P`. Except in the case where
    the GCD is zero, the GCD `G` is made monic.

    Polynomials ``S`` and ``T`` are computed such that
    ``S*A + T*B = G``. The length of ``S`` will be at most
    ``lenB`` and the length of ``T`` will be at most ``lenA``.

.. function:: slong _fq_poly_xgcd(fq_struct *G, fq_struct *S, fq_struct *T, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fmpz_t invB, const fq_ctx_t ctx)

    Computes the GCD of `A` and `B` together with cofactors `S` and `T`
    such that `S A + T B = G`.  Returns the length of `G`.

    Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) \geq 1` and
    `(\operatorname{len}(A),\operatorname{len}(B)) \neq (1,1)`.

    No attempt is made to make the GCD monic.

    Requires that `G` have space for `\operatorname{len}(B)` coefficients.  Writes
    `\operatorname{len}(B)-1` and `\operatorname{len}(A)-1` coefficients to `S` and `T`, respectively.
    Note that, in fact, `\operatorname{len}(S) \leq \max(\operatorname{len}(B) - \operatorname{len}(G), 1)` and
    `\operatorname{len}(T) \leq \max(\operatorname{len}(A) - \operatorname{len}(G), 1)`.

    No aliasing of input and output operands is permitted.

.. function:: void fq_poly_xgcd(fq_poly_t G, fq_poly_t S, fq_poly_t T, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Computes the GCD of `A` and `B`. The GCD of zero polynomials is
    defined to be zero, whereas the GCD of the zero polynomial and some other
    polynomial `P` is defined to be `P`. Except in the case where
    the GCD is zero, the GCD `G` is made monic.

    Polynomials ``S`` and ``T`` are computed such that
    ``S*A + T*B = G``. The length of ``S`` will be at most
    ``lenB`` and the length of ``T`` will be at most ``lenA``.

.. function:: slong _fq_poly_xgcd_euclidean_f(fq_t f, fq_struct *G, fq_struct *S, fq_struct *T, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fmpz_t invB, const fq_ctx_t ctx)

    Either sets `f = 1` and computes the GCD of `A` and `B` together
    with cofactors `S` and `T` such that `S A + T B = G`; otherwise,
    sets `f` to a non-trivial factor of the modulus of ``ctx`` and
    leaves `G`, `S`, and `T` undefined.  Returns the length of `G`.

    Assumes that `\operatorname{len}(A) \geq \operatorname{len}(B) \geq 1` and
    `(\operatorname{len}(A),\operatorname{len}(B)) \neq (1,1)`.

    No attempt is made to make the GCD monic.

    Requires that `G` have space for `\operatorname{len}(B)` coefficients.  Writes
    `\operatorname{len}(B)-1` and `\operatorname{len}(A)-1` coefficients to `S` and `T`, respectively.
    Note that, in fact, `\operatorname{len}(S) \leq \max(\operatorname{len}(B) - \operatorname{len}(G), 1)` and
    `\operatorname{len}(T) \leq \max(\operatorname{len}(A) - \operatorname{len}(G), 1)`.

    No aliasing of input and output operands is permitted.

.. function:: void fq_poly_xgcd_euclidean_f(fq_t f, fq_poly_t G, fq_poly_t S, fq_poly_t T, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

    Either sets `f = 1` and computes the GCD of `A` and `B` or sets
    `f` to a non-trivial factor of the modulus of ``ctx``.

    If the GCD is computed, polynomials ``S`` and ``T`` are
    computed such that ``S*A + T*B = G``; otherwise, they are
    undefined.  The length of ``S`` will be at most ``lenB`` and
    the length of ``T`` will be at most ``lenA``.

    The GCD of zero polynomials is defined to be zero, whereas the GCD
    of the zero polynomial and some other polynomial `P` is defined to
    be `P`. Except in the case where the GCD is zero, the GCD `G` is
    made monic.


Divisibility testing
--------------------------------------------------------------------------------


.. function:: int _fq_poly_divides(fq_struct *Q, const fq_struct *A, slong lenA, const fq_struct *B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

    Returns `1` if ``(B, lenB)`` divides ``(A, lenA)`` exactly and
    sets `Q` to the quotient, otherwise returns `0`.

    It is assumed that `\operatorname{len}(A) \geq \operatorname{len}(B) > 0` and that `Q` has space
    for `\operatorname{len}(A) - \operatorname{len}(B) + 1` coefficients.

    Aliasing of `Q` with either of the inputs is not permitted.

    This function is currently unoptimised and provided for convenience
    only.

.. function:: int fq_poly_divides(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)


    Returns `1` if `B` divides `A` exactly and sets `Q` to the quotient,
    otherwise returns `0`.

    This function is currently unoptimised and provided for convenience
    only.


Derivative
--------------------------------------------------------------------------------


.. function:: void _fq_poly_derivative(fq_struct *rop, const fq_struct *op, slong len, const fq_ctx_t ctx)

    Sets ``(rop, len - 1)`` to the derivative of ``(op, len)``.
    Also handles the cases where ``len`` is `0` or `1` correctly.
    Supports aliasing of ``rop`` and ``op``.

.. function:: void fq_poly_derivative(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

    Sets ``rop`` to the derivative of ``op``.


Evaluation
--------------------------------------------------------------------------------


.. function:: void _fq_poly_evaluate_fq(fq_t rop, const fq_struct *op, slong len, const fq_t a, const fq_ctx_t ctx)

    Sets ``rop`` to ``(op, len)`` evaluated at `a`.

    Supports zero padding.  There are no restrictions on ``len``, that
    is, ``len`` is allowed to be zero, too.

.. function:: void fq_poly_evaluate_fq(fq_t rop, const fq_poly_t f, const fq_t a, const fq_ctx_t ctx)

    Sets ``rop`` to the value of `f(a)`.

    As the coefficient ring `\mathbf{F}_q` is finite, Horner's method
    is sufficient.


Composition
--------------------------------------------------------------------------------


.. function:: void _fq_poly_compose_divconquer(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Computes the composition of ``(op1, len1)`` and ``(op2, len2)``
    using a divide and conquer approach and places the result into ``rop``,
    assuming ``rop`` can hold the output of length
    ``(len1 - 1) * (len2 - 1) + 1``.

    Assumes ``len1, len2 > 0``.  Does not support aliasing between
    ``rop`` and any of ``(op1, len1)`` and ``(op2, len2)``.

.. function:: void fq_poly_compose_divconquer(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the composition of ``op1`` and ``op2``.
    To be precise about the order of composition, denoting ``rop``,
    ``op1``, and ``op2`` by `f`, `g`, and `h`, respectively,
    sets `f(t) = g(h(t))`.

.. function:: void _fq_poly_compose_horner(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``rop`` to the composition of ``(op1, len1)`` and
    ``(op2, len2)``.

    Assumes that ``rop`` has space for ``(len1-1)*(len2-1) + 1``
    coefficients.  Assumes that ``op1`` and ``op2`` are non-zero
    polynomials.  Does not support aliasing between any of the inputs and
    the output.

.. function:: void fq_poly_compose_horner(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the composition of ``op1`` and ``op2``.
    To be more precise, denoting ``rop``, ``op1``, and ``op2``
    by `f`, `g`, and `h`, sets `f(t) = g(h(t))`.

    This implementation uses Horner's method.

.. function:: void _fq_poly_compose(fq_struct *rop, const fq_struct *op1, slong len1, const fq_struct *op2, slong len2, const fq_ctx_t ctx)

    Sets ``rop`` to the composition of ``(op1, len1)`` and
    ``(op2, len2)``.

    Assumes that ``rop`` has space for ``(len1-1)*(len2-1) + 1``
    coefficients.  Assumes that ``op1`` and ``op2`` are non-zero
    polynomials.  Does not support aliasing between any of the inputs and
    the output.

.. function:: void fq_poly_compose(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

    Sets ``rop`` to the composition of ``op1`` and ``op2``.
    To be precise about the order of composition, denoting ``rop``,
    ``op1``, and ``op2`` by `f`, `g`, and `h`, respectively,
    sets `f(t) = g(h(t))`.

.. function:: void _fq_poly_compose_mod_horner(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_ctx_t ctx)


    Sets ``res`` to the composition `f(g)` modulo `h`. We require that
    `h` is nonzero and that the length of `g` is one less than the
    length of `h` (possibly with zero padding). The output is not allowed
    to be aliased with any of the inputs.

    The algorithm used is Horner's rule.

.. function:: void fq_poly_compose_mod_horner(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require that
    `h` is nonzero. The algorithm used is Horner's rule.

.. function:: void _fq_poly_compose_mod_horner_preinv(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_struct * hinv, slong lenhiv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that the length of `g` is one less than
    the length of `h` (possibly with zero padding). We also require
    that the length of `f` is less than the length of
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.  The output is not allowed to be aliased with
    any of the inputs.

    The algorithm used is Horner's rule.

.. function:: void fq_poly_compose_mod_horner_preinv(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_poly_t hinv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that `f` has smaller degree than
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.  The algorithm used is Horner's rule.


.. function:: void _fq_poly_compose_mod_brent_kung(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that the length of `g` is one less than
    the length of `h` (possibly with zero padding). We also require
    that the length of `f` is less than the length of `h`. The output
    is not allowed to be aliased with any of the inputs.

    The algorithm used is the Brent-Kung matrix algorithm.

.. function:: void fq_poly_compose_mod_brent_kung(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that `f` has smaller degree than `h`.  The
    algorithm used is the Brent-Kung matrix algorithm.

.. function:: void _fq_poly_compose_mod_brent_kung_preinv(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_struct * hinv, slong lenhiv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that the length of `g` is one less than
    the length of `h` (possibly with zero padding). We also require
    that the length of `f` is less than the length of
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.  The output is not allowed to be aliased with
    any of the inputs.

    The algorithm used is the Brent-Kung matrix algorithm.

.. function:: void fq_poly_compose_mod_brent_kung_preinv(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_poly_t hinv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that `f` has smaller degree than
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.  The algorithm used is the Brent-Kung matrix
    algorithm.

.. function:: void _fq_poly_compose_mod(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that the length of `g` is one less than
    the length of `h` (possibly with zero padding). The output is not
    allowed to be aliased with any of the inputs.

.. function:: void fq_poly_compose_mod(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero.

.. function:: void _fq_poly_compose_mod_preinv(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_struct * hinv, slong lenhiv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that the length of `g` is one less than
    the length of `h` (possibly with zero padding). We also require
    that the length of `f` is less than the length of
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.  The output is not allowed to be aliased with
    any of the inputs.

.. function:: void fq_poly_compose_mod_preinv(fq_poly_t res, const fq_poly_t f, const fq_poly_t g, const fq_poly_t h, const fq_poly_t hinv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero and that `f` has smaller degree than
    `h`. Furthermore, we require ``hinv`` to be the inverse of the
    reverse of ``h``.

.. function:: void _fq_poly_reduce_matrix_mod_poly (fq_mat_t A, const fq_mat_t B, const fq_poly_t f, const fq_ctx_t ctx)

    Sets the ith row of ``A`` to the reduction of the ith row of `B` modulo
    `f` for `i=1,\ldots,\sqrt{\deg(f)}`. We require `B` to be at least
    a `\sqrt{\deg(f)}\times \deg(f)` matrix and `f` to be nonzero.

.. function:: void _fq_poly_precompute_matrix (fq_mat_t A, const fq_struct* f, const fq_struct* g, slong leng, const fq_struct* ginv, slong lenginv, const fq_ctx_t ctx)

    Sets the ith row of ``A`` to `f^i` modulo `g` for
    `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    be the inverse of the reverse of ``g`` and `g` to be nonzero.

.. function:: void fq_poly_precompute_matrix (fq_mat_t A, const fq_poly_t f, const fq_poly_t g, const fq_poly_t ginv, const fq_ctx_t ctx)

    Sets the ith row of ``A`` to `f^i` modulo `g` for
    `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    be the inverse of the reverse of ``g``.


.. function:: void _fq_poly_compose_mod_brent_kung_precomp_preinv(fq_struct* res, const fq_struct* f, slong lenf, const fq_mat_t A, const fq_struct* h, slong lenh, const fq_struct* hinv, slong lenhinv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero. We require that the ith row of `A` contains
    `g^i` for `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a
    `\sqrt{\deg(h)}\times \deg(h)` matrix. We also require that the
    length of `f` is less than the length of `h`. Furthermore, we
    require ``hinv`` to be the inverse of the reverse of ``h``.
    The output is not allowed to be aliased with any of the inputs.

    The algorithm used is the Brent-Kung matrix algorithm.

.. function:: void fq_poly_compose_mod_brent_kung_precomp_preinv(fq_poly_t res, const fq_poly_t f, const fq_mat_t A, const fq_poly_t h, const fq_poly_t hinv, const fq_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that the ith row of `A` contains `g^i` for
    `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a `\sqrt{\deg(h)}\times
    \deg(h)` matrix. We require that `h` is nonzero and that `f` has
    smaller degree than `h`. Furthermore, we require ``hinv`` to be
    the inverse of the reverse of ``h``. This version of Brent-Kung
    modular composition is particularly useful if one has to perform
    several modular composition of the form `f(g)` modulo `h` for
    fixed `g` and `h`.



Output
--------------------------------------------------------------------------------


.. function:: int _fq_poly_fprint_pretty(FILE *file, const fq_struct *poly, slong len, const char *x, const fq_ctx_t ctx)

    Prints the pretty representation of ``(poly, len)`` to the stream
    ``file``, using the string ``x`` to represent the indeterminate.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int fq_poly_fprint_pretty(FILE * file, const fq_poly_t poly, const char *x, const fq_ctx_t ctx)

    Prints the pretty representation of ``poly`` to the stream
    ``file``, using the string ``x`` to represent the indeterminate.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


.. function:: int _fq_poly_print_pretty(const fq_struct *poly, slong len, const char *x, const fq_ctx_t ctx)

    Prints the pretty representation of ``(poly, len)`` to ``stdout``,
    using the string ``x`` to represent the indeterminate.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


.. function:: int fq_poly_print_pretty(const fq_poly_t poly, const char *x, const fq_ctx_t ctx)

    Prints the pretty representation of ``poly`` to ``stdout``,
    using the string ``x`` to represent the indeterminate.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int _fq_poly_fprint(FILE *file, const fq_struct *poly, slong len, const fq_ctx_t ctx)

    Prints the pretty representation of ``(poly, len)`` to the stream
    ``file``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int fq_poly_fprint(FILE * file, const fq_poly_t poly, const fq_ctx_t ctx)

    Prints the pretty representation of ``poly`` to the stream
    ``file``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


.. function:: int _fq_poly_print(const fq_struct *poly, slong len, const fq_ctx_t ctx)

    Prints the pretty representation of ``(poly, len)`` to ``stdout``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.


.. function:: int fq_poly_print(const fq_poly_t poly, const fq_ctx_t ctx)

    Prints the representation of ``poly`` to ``stdout``.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: char * _fq_poly_get_str(const fq_struct * poly, slong len, const fq_ctx_t ctx)

    Returns the plain FLINT string representation of the polynomial
    ``(poly, len)``.

.. function:: char * fq_poly_get_str(const fq_poly_t poly, const fq_ctx_t ctx)

    Returns the plain FLINT string representation of the polynomial
    ``poly``.

.. function:: char * _fq_poly_get_str_pretty(const fq_struct * poly, slong len, const char * x, const fq_ctx_t ctx)

    Returns a pretty representation of the polynomial
    ``(poly, len)`` using the null-terminated string ``x`` as the
    variable name.

.. function:: char * fq_poly_get_str_pretty(const fq_poly_t poly, const char * x, const fq_ctx_t ctx)

    Returns a pretty representation of the polynomial ``poly`` using the
    null-terminated string ``x`` as the variable name


Inflation and deflation
--------------------------------------------------------------------------------


.. function:: void fq_poly_inflate(fq_poly_t result, const fq_poly_t input, ulong inflation, const fq_ctx_t ctx)

    Sets ``result`` to the inflated polynomial `p(x^n)` where
    `p` is given by ``input`` and `n` is given by ``inflation``.

.. function:: void fq_poly_deflate(fq_poly_t result, const fq_poly_t input, ulong deflation, const fq_ctx_t ctx)

    Sets ``result`` to the deflated polynomial `p(x^{1/n})` where
    `p` is given by ``input`` and `n` is given by ``deflation``.
    Requires `n > 0`.

.. function:: ulong fq_poly_deflation(const fq_poly_t input, const fq_ctx_t ctx)

    Returns the largest integer by which ``input`` can be deflated.
    As special cases, returns 0 if ``input`` is the zero polynomial
    and 1 of ``input`` is a constant polynomial.
