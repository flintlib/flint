.. _gr-poly:

**gr_poly.h** -- dense univariate polynomials over generic rings
===============================================================================

A :type:`gr_poly_t` represents a univariate polynomial `f \in R[X]`
implemented as a dense array of coefficients in a generic ring *R*.

In this module, the context object ``ctx`` always represents the
coefficient ring *R* unless otherwise stated.
Creating a context object representing the polynomial ring `R[X]`
only becomes necessary when one
wants to manipulate polynomials using generic ring methods
like ``gr_add`` instead of the designated polynomial
methods like ``gr_poly_add``.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (often requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Supported coefficient domains
-------------------------------------------------------------------------------

Some methods in this module implicitly assume that *R* is a commutative
ring or an approximate (e.g. floating-point) commutative ring.
When used with a more general *R*, they may output nonsense without
returning the appropriate ``GR_DOMAIN`` or ``GR_UNABLE`` flags.
Better support for noncommutative coefficients is planned for the future.

Some methods make stronger implicit assumptions, for example that *R*
is an integral domain or a field. Such assumptions are documented on
a case by case basis.

Type compatibility
-------------------------------------------------------------------------------

The ``gr_poly`` type has the same data layout as the following
polynomial types: ``fmpz_poly``, ``fq_poly``, ``fq_nmod_poly``,
``fq_zech_poly``, ``arb_poly``, ``acb_poly``, ``ca_poly``.
Methods in this module can therefore be mixed freely with
methods in the corresponding FLINT modules
when the underlying coefficient type is the same.
It is not directly compatible with the following types:
``fmpq_poly`` (coefficients are stored with a common denominator),
``nmod_poly`` (modulus data is stored as part of the polynomial object).

Weak normalization
-------------------------------------------------------------------------------

A :type:`gr_poly_t` is always normalised by removing leading zeros.
For rings without decidable equality (e.g. rings with inexact
representation), only coefficients that are provably zero will be
removed, and there can thus be spurious leading zeros in the
internal representation.
Methods that depend on knowing the exact degree of a polynomial
will act appropriately, typically by returning ``GR_UNABLE``
when it is unknown whether the leading stored coefficient is nonzero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_poly_struct

.. type:: gr_poly_t

    Contains a pointer to an array of coefficients (``coeffs``), the used
    length (``length``), and the allocated size of the array (``alloc``).

    A ``gr_poly_t`` is defined as an array of length one of type
    ``gr_poly_struct``, permitting a ``gr_poly_t`` to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)

.. function:: void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)

.. function:: void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)

.. function:: gr_ptr gr_poly_coeff_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
              gr_srcptr gr_poly_coeff_srcptr(const gr_poly_t poly, slong i, gr_ctx_t ctx)

.. function:: gr_ptr gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
              gr_srcptr gr_poly_entry_srcptr(const gr_poly_t poly, slong i, gr_ctx_t ctx)

    These functions are deprecated aliases of :func:`gr_poly_coeff_ptr` and
    :func:`gr_poly_coeff_srcptr`; use those functions instead.

.. function:: slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: void gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)

.. function:: void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

.. function:: void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
              int gr_poly_get_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx)
              int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx)
              int gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)

.. function:: int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
              int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

.. function:: int gr_poly_truncate(gr_poly_t res, const gr_poly_t poly, slong newlen, gr_ctx_t ctx)

.. function:: int gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_gen(gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_write(gr_stream_t out, const gr_poly_t poly, const char * x, gr_ctx_t ctx)
              int _gr_poly_write(gr_stream_t out, gr_srcptr poly, slong n, const char * x, gr_ctx_t ctx)
              int _gr_poly_get_str(char ** res, const gr_poly_t f, const char * x, gr_ctx_t ctx)
              int gr_poly_get_str(char ** res, const gr_poly_t f, const char * x, gr_ctx_t ctx)
              int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int _gr_poly_set_str(gr_ptr res, const char * s, const char * x, slong len, gr_ctx_t ctx)
              int gr_poly_set_str(gr_poly_t res, const char * s, const char * x, gr_ctx_t ctx)

    Parse polynomial from an expression string, assuming that the string in *x* gives
    the name of the generator. The underscore method zero-pads the result if
    the length of the parsed polynomial is shorter than *len*, and returns
    ``GR_UNABLE`` if the length of the parsed polynomial exceeds *len*.
    Intermediate terms are allowed to be longer than *len*.

    Warning: these methods are not currently optimized for polynomials of high degree
    and may run with quadratic complexity.

.. function:: int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)

.. function:: truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: truth_t gr_poly_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
              truth_t gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)
              truth_t gr_poly_is_gen(const gr_poly_t poly, gr_ctx_t ctx)
              truth_t gr_poly_is_scalar(const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_set_si(gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_set_ui(gr_poly_t poly, ulong c, gr_ctx_t ctx)
              int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
              int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

.. function:: int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong c, gr_ctx_t ctx)
              int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong c, gr_ctx_t ctx)
              int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t c, gr_ctx_t ctx)
              int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t c, gr_ctx_t ctx)

.. function:: int gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

.. function:: int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

.. function:: int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
              int _gr_poly_mullow(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
              int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong len, gr_ctx_t ctx)

    Note: multiplication and low multiplication call :func:`_gr_poly_mullow_generic`
    by default, which currently always delegates to :func:`_gr_poly_mullow_classical`.
    This can be overridden by specific rings.

Multiplication algorithms
-------------------------------------------------------------------------------

.. function:: int _gr_poly_mullow_classical(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_mullow_classical(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

    Multiply using the classical (schoolbook) algorithm, performing
    a sequence of dot products.

.. function:: int _gr_poly_mullow_bivariate_KS(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_mullow_bivariate_KS(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

    Assuming that the coefficients are polynomials of type ``gr_poly``,
    reduce the bivariate product to a single univariate product using
    Kronecker substitution. Returns ``GR_UNABLE`` if the elements are not
    of type ``gr_poly``.

.. function:: int _gr_poly_mullow_complex_reorder(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)
              int gr_poly_mullow_complex_reorder(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, int karatsuba, gr_ctx_t ctx, gr_ctx_t real_ctx)

    Assuming that the coefficients of the polynomials of type ``ctx`` are
    complex numbers represented in Cartesian form as pairs of elements
    of type ``real_ctx``, multiply by reordering to obtain a set of
    real polynomial products.

    If at least one polynomial is detected to be purely real or imaginary, one or two
    real multiplications are used. Otherwise, four real multiplications
    are used unless the *karatsuba* flag is set in which three
    multiplications are used.
    When squaring, two real squarings and one real multiplication are used
    unless the *karatsuba* flag is set in which case three real squarings are used.

.. function:: int _gr_poly_mul_karatsuba(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_mul_karatsuba(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    Karatsuba multiplication.
    Not optimized for unbalanced operands, and not memory-optimized for recursive calls.
    The underscore method requires positive lengths and does not support aliasing.
    This function calls :func:`_gr_poly_mul` recursively rather than itself, so to get a recursive
    algorithm with `O(n^{1.6})` complexity, the ring must overload :func:`_gr_poly_mul` to dispatch
    to :func:`_gr_poly_mul_karatsuba` above some cutoff.

.. function:: int _gr_poly_mul_toom33(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
              int gr_poly_mul_toom33(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

    Balanced Toom-3 multiplication with interpolation in five points,
    using the Bodrato evaluation scheme. Assumes commutativity and that the ring
    supports exact division by 2 and 3.
    Not optimized for squaring.
    The underscore method requires positive lengths and does not support aliasing.
    This function calls :func:`_gr_poly_mul` recursively rather than itself, so to get a recursive
    algorithm with `O(n^{1.465})` complexity, the ring must overload :func:`_gr_poly_mul` to dispatch
    to :func:`_gr_poly_mul_toom33` above some cutoff.

.. function:: int _gr_poly_mullow_toom_serial(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong trunc, slong num_points, slong splitting_degree, gr_ctx_t ctx)
              int gr_poly_mullow_toom_serial(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, slong num_points, slong splitting_degree, gr_ctx_t ctx)

    Toom-Cook multiplication of arbitrary order: the inputs are segmented as
    polynomials in `t = x^m` for some splitting degree `m \ge 1` and the product
    is computed by evaluation and interpolation at `r \ge 1` distinct
    values of `t`.

    For convenience, the user can specify either the number of interpolation
    points `r` (``num_points``) or the splitting degree  `m` (``splitting_degree``).
    The user must provide exactly one of these parameters, passing zero for the
    other; the optimal value of the omitted parameter will be computed automatically.
    For example, if ``len1`` and ``len2`` are roughly balanced, then providing
    either `r = 3, m = 0` or `r = 0, m = \lceil \max(len1, len2) / 2 \rceil`
    will result in Karatsuba multiplication with evaluation in `0`, `1` and `\infty`.
    For balanced operands, an odd number of evaluation points `r = 2k-1`
    corresponds to splitting both inputs into `k` segments; even values of `d`
    typically only make sense for unbalanced products.

    The evaluation and interpolation points are chosen as
    `0, 1, -1, 2, -2, 3, -3, \ldots, \infty`. Evaluation and interpolation is
    implicitly done by multiplying by the corresponding Vandermonde matrix
    and its inverse. In too small finite characteristic where the Vandermonde
    matrix cannot be inverted, ``GR_UNABLE`` is returned. In characteristic zero,
    the interpolation matrix is computed with a common denominator and this
    denominator is removed from the final result with an exact division.


    The suffix "serial" indicates that the evaluations are done one point at a time
    and that the output is constructed incrementally with the goal to minimize
    memory allocation. This implementation is therefore well suited for huge
    memory-constrained multiplications, in contrast to a typical Toom-Cook
    implementation which will use several extra temporary buffers to
    minimize scalar operations. In total this implementation requires `4m`
    elements of temporary storage for intermediate polynomials and `2r^2`
    ring elements or integers for the evaluation and interpolation matrices.

    The underscore method requires positive lengths and does not support aliasing.

.. function:: int gr_poly_add_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_add_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
              int gr_poly_add_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_add_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz c, gr_ctx_t ctx)
              int gr_poly_add_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq c, gr_ctx_t ctx)

    Sets *res* to *poly* plus the scalar *c* which must be
    an element of or coercible to the coefficient ring.

.. function:: int gr_poly_sub_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_sub_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
              int gr_poly_sub_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_sub_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz c, gr_ctx_t ctx)
              int gr_poly_sub_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq c, gr_ctx_t ctx)

    Sets *res* to *poly* minus the scalar *c* which must be
    an element of or coercible to the coefficient ring.

.. function:: int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_scalar_mul(gr_poly_t res, gr_srcptr c, const gr_poly_t poly, gr_ctx_t ctx)
              int gr_poly_mul_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
              int gr_poly_mul_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
              int gr_poly_mul_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz c, gr_ctx_t ctx)
              int gr_poly_mul_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq c, gr_ctx_t ctx)

    Sets *res* to *poly* multiplied by the scalar *c* (or the scalar *c* multiplied by *poly*)
    which must be an element of or coercible to the coefficient ring.

.. function:: int gr_poly_addmul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

    Adds *poly* multiplied by the scalar *c* to *res*.

.. function:: int gr_poly_submul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

    Subtracts *poly* multiplied by the scalar *c* from *res*.

Powering
--------------------------------------------------------------------------------

.. function:: int _gr_poly_pow_series_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
              int gr_poly_pow_series_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_series_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
              int gr_poly_pow_series_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
              int gr_poly_pow_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
              int gr_poly_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)

.. function:: int gr_poly_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx)

.. function:: int _gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t exp, slong len, int precomp, gr_ctx_t ctx)
              int gr_poly_pow_series_fmpq_recurrence(gr_poly_t res, const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx)

Shifting
-------------------------------------------------------------------------------

.. function:: int _gr_poly_shift_left(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
              int gr_poly_shift_left(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

.. function:: int _gr_poly_shift_right(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
              int gr_poly_shift_right(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)


Scalar division
--------------------------------------------------------------------------------

.. function:: int gr_poly_div_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

Division with remainder
--------------------------------------------------------------------------------

.. function:: int _gr_poly_divrem_divconquer_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_divrem_divconquer_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
              int gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int _gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    These functions implement Euclidean division with remainder:
    given polynomials `A, B \in K[x]` where `K` is a field, with `B \ne 0`,
    there is a unique quotient `Q` and remainder `R` such that `A = BQ + R`
    and either `R = 0` or `\deg(R) < \deg(B)`.
    If *B* is provably zero, ``GR_DOMAIN`` is returned.

    When `K` is a commutative ring and `\operatorname{lc}(B)` is a unit in `K`,
    the situation is the same as over fields. In particular, Euclidean division
    with remainder always makes sense over commutative rings when `B` is monic.
    If `\operatorname{lc}(B)` is not a unit, the division still makes sense if
    the coefficient quotient `\operatorname{lc}(r)`  / `\operatorname{lc}(B)`
    exists for each partial remainder `r`. Indeed,
    the *basecase* and *divconquer* algorithms return ``GR_DOMAIN`` precisely when
    encountering a leading quotient `\operatorname{lc}(r)`  / `\operatorname{lc}(B) \not \in K`.
    However, the *newton* algorithm as currently implemented
    returns ``GR_DOMAIN`` when `\operatorname{lc}(B)^{-1} \not \in K`.

    The underscore methods make the following assumptions:

    * *Q* has room for ``lenA - lenB + 1`` coefficients.
    * *R* has room for ``lenB - 1`` coefficients.
    * ``lenA >= lenB >= 1``.
    * *Q* is not aliased with either *A* or *B*.
    * *R* is not aliased with *B*.
    * *R* may be aliased with *A*, in which case all ``lenA``
      entries may be used as scratch space. Note that in this case,
      only the low ``lenB - 1`` coefficients of *R* actually represent
      valid coefficients on output: the higher scratch coefficients will not
      necessarily be zeroed.
    * The divisor *B* is normalized to have nonzero leading coefficient.
      (The non-underscore methods check for leading coefficients that
      are not provably nonzero and return ``GR_UNABLE``.)

    The *preinv1* functions take a precomputed inverse of the
    leading coefficient as input.
    The *noinv* versions perform repeated checked divisions
    by the leading coefficient.

.. function:: int _gr_poly_div_divconquer_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_div_divconquer_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_div_divconquer(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
              int gr_poly_div_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_div_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
              int _gr_poly_div_basecase_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int _gr_poly_div_basecase(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_div_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_div_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_div(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_div(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    Versions of the *divrem* functions which output only the quotient.
    These are generally faster.

.. function:: int _gr_poly_rem(gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    Versions of the *divrem* functions which output only the remainder.

Division with remainder with full precomputed inverse
--------------------------------------------------------------------------------

.. function:: int _gr_poly_div_newton_n_preinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr FLINT_UNUSED(B), slong lenB, gr_srcptr Binv, slong lenBinv, gr_ctx_t ctx)

    Notionally computes polynomials `Q` and `R` such that `A = BQ + R` with
    `\operatorname{len}(R)` less than ``lenB``, where ``A`` is of length ``lenA``
    and ``B`` is of length ``lenB``, but return only `Q`.

    We require that `Q` have space for ``lenA - lenB + 1`` coefficients
    and assume that the leading coefficient of `B` is a unit. Furthermore, we
    assume that `Binv` is the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.

    The algorithm used is to reverse the polynomials and divide the
    resulting power series, then reverse the result.

.. function:: int gr_poly_div_newton_n_preinv(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, const gr_poly_t Binv, gr_ctx_t ctx)

    Notionally computes `Q` and `R` such that `A = BQ + R` with
    `\operatorname{len}(R) < \operatorname{len}(B)`, but returns only `Q`.

    We assume that the leading coefficient of `B` is a unit and that `Binv` is
    the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`.

    It is required that the length of `A` is less than or equal to
    2*the length of `B` - 2.

    The algorithm used is to reverse the polynomials and divide the
    resulting power series, then reverse the result.

.. function:: int _gr_poly_divrem_newton_n_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr Binv, slong lenBinv, gr_ctx_t ctx)

    Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R)` less
    than ``lenB``, where `A` is of length ``lenA`` and `B` is of
    length ``lenB``. We require that `Q` have space for
    ``lenA - lenB + 1`` coefficients. Furthermore, we assume that `Binv` is
    the inverse of the reverse of `B` mod `x^{\operatorname{len}(B)}`. The algorithm
    used is to call :func:`_gr_poly_div_newton_n_preinv` and then multiply out
    and compute the remainder.

.. function:: int gr_poly_divrem_newton_n_preinv(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, const gr_poly_t Binv, gr_ctx_t ctx)

    Computes `Q` and `R` such that `A = BQ + R` with `\operatorname{len}(R) <
    \operatorname{len}(B)`.  We assume `Binv` is the inverse of the reverse of `B`
    mod `x^{\operatorname{len}(B)}`.

    It is required (not checked) that the length of `A` is less than or equal to
    2*the length of `B` - 2.


Power series division
--------------------------------------------------------------------------------

For divide-and-conquer (including Newton-like) algorithms, *cutoff* has the
following meaning: we use the basecase algorithm for lengths `n < \operatorname{cutoff}`
and the divide-and-conquer algorithm for `n \ge \operatorname{cutoff}`.
Using `\operatorname{cutoff} = n` thus results in exactly one divide-and-conquer
step with a basecase length of `\lceil n / 2 \rceil`.
One should **avoid** calling the Newton methods with `n < \operatorname{cutoff}`
as this may result in much worse performance if those methods
do not have a specific escape check for that case.

The *newton* versions uses Newton iteration, switching to a basecase
algorithm when the length is smaller than the specified *cutoff*.
Division uses the Karp-Markstein algorithm.

.. function:: int _gr_poly_inv_series_newton(gr_ptr res, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_inv_series_newton(gr_poly_t res, const gr_poly_t A, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_inv_series_basecase_preinv1(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr Ainv, slong len, gr_ctx_t ctx)
              int _gr_poly_inv_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_inv_series_basecase(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
              int _gr_poly_inv_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_inv_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_div_series_newton(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_div_series_newton(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_div_series_divconquer(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_div_series_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_div_series_invmul(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series_invmul(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_srcptr Binv, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series_basecase(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
              int _gr_poly_div_series(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_div_series(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)

Exact division
--------------------------------------------------------------------------------

These functions compute a quotient `Q = A / B` which is known to be exact
(without remainder) in `R[x]` (or in `R[[x]] / x^n` in the case of series
division). Given a nonexact division, they are allowed to set `Q` to
an arbitrary polynomial and return ``GR_SUCCESS`` instead of returning an
error flag.

`R` is assumed to be an integral domain (this is not checked).

For exact division, we have the choice of starting the division
from the most significant terms (classical division) or the least significant
(power series division). Which direction is more efficient depends
in part on whether the leading or trailing coefficient of `B` is cheaper
to use for divisions. In a generic setting, this is hard to predict.

The *bidirectional* algorithms combine two half-divisions from both ends.
This halves the number of operations in the basecase regime, though an
extra coefficient inversion may be needed.

The ``noinv`` versions perform repeated ``divexact`` operations in the
scalar domain without attempting to invert the leading (or trailing) coefficient,
while other versions check invertibility first.
There are no ``divexact_preinv1`` versions because those are identical to the
``div_preinv1`` counterparts.

.. function:: int _gr_poly_divexact_basecase_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
              int gr_poly_divexact_basecase_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_divexact_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
              int gr_poly_divexact_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_divexact_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
              int _gr_poly_divexact_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
              int gr_poly_divexact_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

.. function:: int _gr_poly_divexact_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int _gr_poly_divexact_series_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
              int gr_poly_divexact_series_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)


Square roots
--------------------------------------------------------------------------------

.. function:: int _gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_sqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_rsqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_rsqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

Evaluation
-------------------------------------------------------------------------------

.. function:: int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_modular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate_modular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

    Set *res* to *poly* evaluated at *x*.

.. function:: int _gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int _gr_poly_evaluate_other_rectangular(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int _gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    Set *res* to *poly* evaluated at *x*, where the coefficients of *f*
    belong to *ctx* while both *x* and *res* belong to *x_ctx*.


Newton basis
-------------------------------------------------------------------------------

Let `C = (c_0, c_1, c_2, \ldots)` be a sequence of
values of the commutative ring *R*.
Then the sequence of polynomials

.. math::

      \begin{matrix}
      N_0 & = & 1, \\
      N_1 & = & (x-c_0), \\
      N_2 & = & (x-c_0)(x-c_1), \\
      N_3 & = & (x-c_0)(x-c_1)(x-c_2), \ldots
      \end{matrix}

with roots in `C` forms a basis of `R[x]` called the *Newton basis*
associated with *C*.
Given a polynomial `f = f_0 + f_1 x + \ldots + f_n x^n` in the
standard monomial basis, there are accordingly coefficients `g_0, \ldots, g_n`
such that `f = g_0 N_0 + g_1 N_1 + \ldots + g_n N_n`
and vice versa.
In the following functions, some finite initial segment of `C` is represented
by a vector called *basis*, and a polynomial in the Newton basis is
represented formally as the polynomial
`g = g_0 + g_1 x + \ldots + g_n x^n`.

The following functions all support aliasing between inputs and outputs
except that no output may be aliased with *basis*.

.. function:: int _gr_poly_newton_basis_from_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_newton_basis_from_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)
              int _gr_poly_newton_basis_to_monomial(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_newton_basis_to_monomial(gr_poly_t res, const gr_vec_t basis, const gr_poly_t poly, gr_ctx_t ctx)

    Given *poly* of length *len* in the standard monomial basis, set *res* to
    the polynomial in the Newton basis and vice versa.
    The underscore methods require that ``basis`` is a pointer to at least
    `len - 1` entries. The non-underscore methods return ``GR_UNABLE`` if
    *basis* does not contain at least `len - 1` entries.

.. function:: int _gr_poly_newton_basis_evaluate(gr_ptr res, gr_srcptr basis, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_poly_newton_basis_evaluate(gr_ptr res, const gr_vec_t basis, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

    Given *poly* of length *len* in the Newton basis, set *res* to the
    evaluation at *x*.
    The underscore method requires that ``basis`` is a pointer to at least
    `len - 1` entries. The non-underscore method returns ``GR_UNABLE`` if
    *basis* does not contain at least `len - 1` entries.

.. function:: int _gr_poly_newton_basis_interpolate_exact(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
              int gr_poly_newton_basis_interpolate_exact(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)
              int _gr_poly_newton_basis_interpolate(gr_ptr res, gr_srcptr basis, gr_srcptr ys, slong len, gr_ctx_t ctx)
              int gr_poly_newton_basis_interpolate(gr_poly_t res, const gr_vec_t basis, const gr_vec_t ys, gr_ctx_t ctx)

    Given a vector *y* containing *len* values, set *res* to the Newton basis form
    of the unique interpolating
    polynomial *f* of length (up to) *len* such that
    `y_0 = f(c_0), \ldots, y_{len-1} = f(c_{len-1})`
    evaluated at the roots of the Newton basis.

    We presume, but do not check, that the ring is an integral domain.
    These functions may succeed over non-integral domains, but
    the result need not be an interpolating polynomial.

    We require that *basis* contains at least *len* initial roots which
    are pairwise distinct.

    These functions return ``GR_DOMAIN`` (and/or ``GR_UNABLE``)
    in either of the following situations:

    * The vector *basis* does not have sufficiently many points (checked by the
      non-underscore functions only).
    * The evaluation points *xs* are not pairwise distinct.
    * The interpolating polynomial *f* has coefficients in the fraction field
      of *R* but not in *R* itself.

    The *exact* versions presume that the evaluation points are distinct
    and that *f* has coefficients in *R*; they may silently output some
    arbitary polynomial otherwise.


Multipoint evaluation and interpolation
-------------------------------------------------------------------------------

.. function:: gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx)

.. function:: void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx)

    Initialize a subproduct tree over the given roots.

.. function:: int _gr_poly_product_roots(gr_ptr poly, gr_srcptr xs, slong n, gr_ctx_t ctx)
              int gr_poly_product_roots(gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

    Set *poly* to the polynomial `(x-x_0) (x-x_1) \cdots (x-x_{n-1})`.

.. function:: int _gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly, slong plen, const gr_ptr * tree, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
              int gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

.. function:: int _gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
              int gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

.. function:: int _gr_poly_interpolate_exact(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
              int gr_poly_interpolate_exact(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
              int _gr_poly_interpolate(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
              int gr_poly_interpolate(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)

    Given vectors *xs* and *ys* of length *len* where the entries of *xs*
    are pairwise distinct, set *res* to the
    interpolating polynomial *f* of length (up to) *len* such that
    `y_0 = f(x_0), \ldots, y_{len-1} = f(x_{len-1})`.

    We presume, but do not check, that the ring is an integral domain.
    These functions may succeed over non-integral domains, but
    the result might not be an interpolating polynomial.

    These functions return ``GR_DOMAIN`` (and/or ``GR_UNABLE``)
    in either of the following situations:

    * The vectors *xs* and *ys* do not have the same length (checked by the
      non-underscore functions only).
    * The evaluation points *xs* are not pairwise distinct.
    * The interpolating polynomial *f* has coefficients in the fraction field
      of *R* but not in *R* itself.

    The *exact* versions presume that the evaluation points are distinct
    and that *f* has coefficients in *R*; they may silently output some
    arbitary polynomial otherwise.

.. function:: int _gr_poly_interpolation_weights(gr_ptr w, const gr_ptr * tree, slong len, gr_ctx_t ctx)
              int _gr_poly_interpolate_fast_precomp(gr_ptr poly, gr_srcptr ys, const gr_ptr * tree, gr_srcptr weights, slong len, gr_ctx_t ctx)
              int _gr_poly_interpolate_fast(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
              int gr_poly_interpolate_fast(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)

    Fast polynomial interpolation using a subproduct tree. The *precomp*
    version requires a precomputed subproduct tree generated using
    :func:`_gr_poly_tree_build` and precomputed interpolation weights
    generated using :func:`_gr_poly_interpolation_weights`.

    This currently requires a field.


Composition
-------------------------------------------------------------------------------

.. function:: int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int _gr_poly_taylor_shift_convolution(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift_convolution(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
              int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

    Sets *res* to the Taylor shift `f(x+c)`, where *f* is given by
    *poly*, computed respectively using
    an optimized form of Horner's rule, divide-and-conquer, a single
    convolution, and an automatic choice between the three algorithms.
    The underscore methods support aliasing.

.. function:: int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
              int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
              int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    Sets *res* to the composition `f(g(x))` where *f* is given by *poly1*
    and *g* is given by *poly2*, respectively using Horner's rule,
    divide-and-conquer, and an automatic choice between the two algorithms.
    The default algorithm also handles special-form input `g = ax^n + c`
    efficiently by performing a Taylor shift followed by a rescaling.
    The underscore methods do not support aliasing of the output
    with either input polynomial.

Power series composition and reversion
-------------------------------------------------------------------------------

.. function:: int _gr_poly_compose_series_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_compose_series_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
              int _gr_poly_compose_series_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_compose_series_brent_kung(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
              int _gr_poly_compose_series_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_compose_series_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
              int _gr_poly_compose_series(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
              int gr_poly_compose_series(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)

    Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    respectively using Horner's rule, the Brent-Kung baby step-giant step
    algorithm [BrentKung1978]_, divide-and-conquer, and an automatic choice between the algorithms.

    The default algorithm also handles short input and
    special-form input `g = ax^n` efficiently.

    We require that the constant term in `g(x)` is exactly zero.
    The underscore methods do not support aliasing of the output
    with either input polynomial, and do not zero-pad the result.

.. function:: int _gr_poly_revert_series_lagrange(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
              int gr_poly_revert_series_lagrange(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
              int _gr_poly_revert_series_lagrange_fast(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
              int gr_poly_revert_series_lagrange_fast(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
              int _gr_poly_revert_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
              int gr_poly_revert_series_newton(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
              int _gr_poly_revert_series(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
              int gr_poly_revert_series(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)

    Sets *res* to the power series reversion `f^{-1}(x)` which satisfies
    `f^{-1}(f(x)) = f(f^{-1}(x)) = x` mod `x^n`.
    For the series reversion to exist, we require that the constant term
    in `f` is zero and that the linear coefficient is invertible.
    The flag ``GR_DOMAIN`` is returned otherwise.

    The *lagrange* and *lagrange_fast* algorithms require the ability
    to divide by `2, 3, \ldots, n-1` and will return
    the ``GR_UNABLE`` flag in too small characteristic.

    The underscore methods do not support aliasing of the output
    with the input.

    The Newton method is described in [BrentKung1978]_; the
    *lagrange* algorithm implements the Lagrange inversion formula,
    while the *lagrange_fast* algorithm implements the baby-step
    giant-step algorithm described in [Joh2015b]_.

Derivative and integral
-------------------------------------------------------------------------------

.. function:: int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

.. function:: int _gr_poly_nth_derivative(gr_ptr res, gr_srcptr poly, ulong n, slong len, gr_ctx_t ctx)
              int gr_poly_nth_derivative(gr_poly_t res, const gr_poly_t poly, ulong n, gr_ctx_t ctx)

.. function:: int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

Monic polynomials
-------------------------------------------------------------------------------

.. function:: int _gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
              int gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

.. function:: truth_t _gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx)
              truth_t gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx)

.. function:: int gr_poly_canonical_associate(gr_poly_t res, gr_poly_t u, const gr_poly_t src, gr_ctx_t ctx)

GCD
-------------------------------------------------------------------------------

.. function:: int _gr_poly_hgcd(gr_ptr r, slong * sgn, gr_ptr * M, slong * lenM, gr_ptr A, slong * lenA, gr_ptr B, slong * lenB, gr_srcptr a, slong lena, gr_srcptr b, slong lenb, slong cutoff, gr_ctx_t ctx)

    Computes the HGCD of `a` and `b`, that is, a matrix `M`, a sign `\sigma`
    and two polynomials `A` and `B` such that

    .. math::

        (A,B)^t = \sigma M^{-1} (a,b)^t.

    Assumes that `\operatorname{len}(a) > \operatorname{len}(b) > 0`.

    Assumes that `A` and `B` have space of size at least `\operatorname{len}(a)`
    and `\operatorname{len}(b)`, respectively.  On exit, ``*lenA`` and ``*lenB``
    will contain the correct lengths of `A` and `B`.

    Assumes that ``M[0]``, ``M[1]``, ``M[2]``, and ``M[3]``
    each point to a vector of size at least `\operatorname{len}(a)`.

    If `r` is not ``NULL``, writes to that variable the corresponding value
    for computing resultants using the HGCD algorithm.

.. function:: int _gr_poly_gcd_hgcd(gr_ptr G, slong * _lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
              int gr_poly_gcd_hgcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_gcd_subresultant(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_gcd_subresultant(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
              int _gr_poly_gcd_generic(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int _gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    GCD in the polynomial ring `R[x]`.

    The *gcd_subresultant* algorithm assumes that *R* is a unique factorization
    domain. The *euclidean* and *hgcd* algorithms assume that *R* is a field.
    The time complexity of the half-GCD algorithm is `\mathcal{O}(n \log^2 n)`
    field operations. For further details, see [ThullYap1990]_.

    The *gcd_generic* fallback implementation checks if *R* is a field
    or UFD and dispatches to an appropriate algorithm, and otherwise
    returns ``GR_UNABLE``. The main *gcd* function is synonymous with
    *gcd_generic* unless overridden by a particular ring.

    The non-underscore methods canonicalise the leading coefficient of the output
    by calling ``gr_poly_canonical_associate``. Over a field, this corresponds
    to making the output monic. The underscore methods need not perform this
    canonicalisation.

    The underscore methods assume ``lenA >= lenB >= 1`` and that both
    *A* and *B* have nonzero leading coefficient.

.. function:: int _gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_xgcd_euclidean(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

.. function:: int _gr_poly_xgcd_hgcd(slong * Glen, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)
              int gr_poly_xgcd_hgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)

.. function:: int _gr_poly_xgcd_generic(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int _gr_poly_xgcd(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
              int gr_poly_xgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

Resultant
-------------------------------------------------------------------------------

For two non-zero polynomials `f(x) = a_m x^m + \dotsb + a_0` and
`g(x) = b_n x^n + \dotsb + b_0` of degrees `m` and `n`, the resultant
is defined to be

.. math::

        a_m^n b_n^m \prod_{(x, y) : f(x) = g(y) = 0} (x - y).

For convenience, we define the resultant to be equal to zero if either
of the two polynomials is zero.

.. function:: int _gr_poly_resultant_euclidean(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_resultant_euclidean(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
              int _gr_poly_resultant_hgcd(gr_ptr res, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
              int gr_poly_resultant_hgcd(gr_ptr res, const gr_poly_t f, const gr_poly_t g, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_resultant_sylvester(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_resultant_sylvester(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
              int _gr_poly_resultant_small(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_resultant_small(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
              int _gr_poly_resultant(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
              int gr_poly_resultant(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)

    Sets *res* to the resultant of *poly1* and *poly2*.
    The underscore methods assume that `len1 \ge len2 \ge 1`
    and that the leading coefficients are nonzero.

    The *euclidean* algorithm is the ordinary Euclidean algorithm.
    The *hgcd* version uses the quasilinear half-GCD algorithm.
    It requires two extra tuning parameters ``inner_cutoff``
    (recursion threshold passed forward to the HGCD algorithm)
    and ``cutoff``. Both algorithms can fail when run over
    non-fields; they will return ``GR_DOMAIN``
    when encountering an impossible inverse.

    The *small* version uses division-free straight-line programs
    optimized for short polynomials.
    It returns ``GR_UNABLE`` if the polynomials are too large.
    Currently this function handles the cases where `len1 \le 2`
    or `len2 \le 3`.

    The *sylvester* version constructs the Sylvester matrix
    and computes its determinant. This is useful over inexact rings
    and as a fallback for rings without division.

    The default version attempts to choose an appropriate
    algorithm automatically.

    Currently no algorithm has been implemented that is appropriate for
    integral domains.


Squarefree factorization
-------------------------------------------------------------------------------

TODO: currently only fields of characteristic 0 are supported.

.. function:: int gr_poly_factor_squarefree(gr_ptr c, gr_vec_t fac, gr_vec_t exp, const gr_poly_t poly, gr_ctx_t ctx)

    Computes a squarefree factorization of *poly*.

    The constant *c* is set to an element of the scalar ring.
    The factors in *fac* are set to polynomials; the user must thus
    initialize it to a vector of polynomials of the same type as
    *poly* (and *not* to the parent *ctx*).
    The exponent vector *exp* must be initialized to the *fmpz* type.

.. function:: int gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

    Sets *res* to the squarefreepart of *poly*.

Shift equivalence
-------------------------------------------------------------------------------

.. function:: truth_t gr_poly_shift_equivalent(fmpz_t shift, const gr_poly_t p, const gr_poly_t q, gr_ctx_t ctx)

    Returns whether there exists an integer *n* such that `p(x + n) = q(x)`. If
    the result is ``T_TRUE`` and *shift* is not ``NULL``, *shift* is set to
    such an *n*.

.. function:: int gr_poly_leading_taylor_shift(gr_ptr shift, const gr_poly_t p, const gr_poly_t q, gr_ctx_t ctx)

    Computes (if possible) *s* such that `p(x+s) = q(x)(1+O(x^2))`.

.. function:: int gr_poly_dispersion_resultant(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);
              int gr_poly_dispersion_factor(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);
              int gr_poly_dispersion(fmpz_t disp, gr_vec_t disp_set, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);

    Computes the dispersion and/or the dispersion set of *f* and *g*.

    The dispersion set of two polynomials *f* and *g* (over a unique
    factorization domain of characteristic zero) is the set of nonnegative
    integers *n* such that `f(x + n)` and `g(x)` have a nonconstant common
    factor. The dispersion is the largest element of the dispersion set.

    The output variables *disp* and/or *disp_set* can be ``NULL``, in which case
    the corresponding result is not stored.
    When the dispersion set is empty, *disp* is left unchanged.
    The elements of *disp_set* are sorted in increasing order.

    The *factor* version uses the algorithm described in [ManWright1994]_.
    The *resultant* version computes the integer roots of a bivariate resultant
    and is mainly intended for testing.

.. function:: int gr_poly_dispersion_from_factors(fmpz_t disp, gr_vec_t disp_set, const gr_vec_t ffac, const gr_vec_t gfac, gr_ctx_t ctx);

    Same as :func:`gr_poly_dispersion_factor` for nonzero *f* and *g* but takes
    as input their nonconstant irreducible factors (without multiplicities)
    instead of the polynomials themselves.

.. function:: int gr_poly_shiftless_decomposition_factor(gr_ptr c, gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_t f, gr_ctx_t ctx)
              int gr_poly_shiftless_decomposition(gr_ptr c, gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_poly_t f, gr_ctx_t ctx)


    Computes a decomposition of *f* of the form

        .. math:: c \prod_i \prod_j g_i(x + h_{i,j})^{e_{i,j}}

    where

    * `c` is a constant,
    * the `g_i` are squarefree polynomials of degree at least one,
    * `g_i(x)` and `g_j(x + h)` (with `i \neq j`) are coprime for all
      `h \in \mathbb Z`,
    * `g_i(x)` and `g_i(x + h)` are coprime for all nonzero `h \in \mathbb Z`,
    * `e_{i,j}` and `h_{i,j}` are integers with `e_{i,j} \geq 1`
      and `0 = h_{i,1} < h_{i,2} < \cdots`.

    The output variable *slfac* must be initialized to a vector of polynomials
    of the same type as *f*. The other two output vectors *slshift* and
    *slmult* must be initialized to vectors *of vectors* with entries of type
    *fmpz*.

    The *factor* version computes an irreducible factorization and sorts the
    factors into shift-equivalence classes.

    No algorithm avoiding a full irreducible factorization is currently
    implemented.

.. function:: int _gr_poly_shiftless_decomposition_from_factors(gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_vec_t fac, const gr_vec_t mult, gr_ctx_t ctx)
              int gr_poly_shiftless_decomposition_from_factors(gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult, const gr_vec_t fac, const gr_vec_t mult, gr_ctx_t ctx)

    Same as :func:`gr_poly_shiftless_decomposition_factor` but takes as input
    an irreducible factorization (*fac*, *mult*) of *f* (without the
    prefactor *c*). The underscore method does not support aliasing of *slfac*
    with *fac*.

Roots
-------------------------------------------------------------------------------

.. function:: int gr_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
              int gr_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t poly_ctx, int flags, gr_ctx_t ctx)

    Finds all roots of the given polynomial in the ring defined by *ctx*,
    storing the roots without duplication in *roots* (a vector with
    elements of type ``ctx``) and the corresponding multiplicities in
    *mult* (a vector with elements of type ``fmpz``).

    If the target ring is not an algebraically closed field, then
    the sum of multiplicities can be smaller than the degree of the
    polynomial. For example, with ``fmpz`` coefficients, we only
    find integer roots.
    The *other* version of this function takes as input a polynomial
    with entries in a different ring ``poly_ctx``. For example,
    we can compute ``qqbar`` or ``arb`` roots for a polynomial
    with ``fmpz`` coefficients.

    Whether the roots are sorted in any particular order is
    ring-dependent.

    We consider roots of the zero polynomial to be ill-defined and return
    ``GR_DOMAIN`` in that case.

.. function:: int _gr_poly_refine_roots_aberth(gr_ptr w, gr_srcptr f, gr_srcptr f_prime, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)
              int _gr_poly_refine_roots_wdk(gr_ptr w, gr_srcptr f, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)

    Given a vector of approximate complex roots `z_1, \ldots, z_{deg}`
    of `f = \sum_{i=0}^{deg} f_i x^i`,
    computes a vector of corrections `w_1, \ldots, w_{deg}` such that
    `z_k - w_k` is a closer approximation of the respective root
    provided that the initial approximations are close enough
    and that the polynomial evaluation is numerically accurate.
    The user will typically call these methods in a loop.

    The *wdk* version performs the Weierstrass-Durand-Kerner update

    .. math ::

        w_k = \frac{f(z_k)}{\prod_{j \ne k} (z_k - z_j)}, \quad k = 1, \ldots, deg.

    The *aberth* version performs the Aberth-Ehrlich update

    .. math ::

        w_k = \frac{g(z_k)}{1 - g(z_k) \sum_{j \ne k} (z_k - z_j)^{-1}}, \quad g(z_k) = \frac{f(z_k)}{f'(z_k)} \quad k = 1, \ldots, deg.

    requiring the coefficients of `f'` as an extra input *f_prime*.

    If *progressive* flag is set, corrected roots `z_j - w_j` that
    have already been computed are used in place of `z_j` in the
    update loop, which can improve the rate of convergence.

Power series special functions
--------------------------------------------------------------------------------

.. function:: int _gr_poly_asin_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_asin_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_asinh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_asinh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_acos_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_acos_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_acosh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_acosh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_atan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_atan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_atanh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_atanh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
              int _gr_poly_log1p_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_log1p_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

.. function:: int _gr_poly_exp_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
              int gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
              int _gr_poly_exp_series_basecase_mul(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
              int gr_poly_exp_series_basecase_mul(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
              int _gr_poly_exp_series_newton(gr_ptr f, gr_ptr g, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx)
              int gr_poly_exp_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_exp_series_generic(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
              int _gr_poly_exp_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
              int gr_poly_exp_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)

.. function:: int _gr_poly_sin_cos_series_basecase(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
              int gr_poly_sin_cos_series_basecase(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
              int _gr_poly_sin_cos_series_tangent(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
              int gr_poly_sin_cos_series_tangent(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)

    The *basecase* version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The *tangent* version uses the tangent half-angle formulas to compute
    the sine and cosine via :func:`_acb_poly_tan_series`. This
    requires `O(M(n))` operations.
    When `h = h_0 + h_1` where the constant term `h_0` is nonzero,
    the evaluation is done as
    `\sin(h_0 + h_1) = \cos(h_0) \sin(h_1) + \sin(h_0) \cos(h_1)`,
    `\cos(h_0 + h_1) = \cos(h_0) \cos(h_1) - \sin(h_0) \sin(h_1)`.

    The *basecase* and *tangent* versions take a flag *times_pi*
    specifying that the input is to be multiplied by `\pi`.

.. function:: int _gr_poly_tan_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
              int gr_poly_tan_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
              int _gr_poly_tan_series_newton(gr_ptr f, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx)
              int gr_poly_tan_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)
              int _gr_poly_tan_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
              int gr_poly_tan_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)

Modular arithmetic and composition
--------------------------------------------------------------------------------

.. function:: int _gr_poly_mulmod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_srcptr f, slong lenf, gr_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

    It is required that ``len1 + len2 - lenf > 0``, which is
    equivalent to requiring that the result will actually be
    reduced. Otherwise, simply use ``_fq_poly_mul`` instead.

    Aliasing of ``f`` and ``res`` is not permitted.

.. function:: int gr_poly_mulmod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t f, gr_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

.. function:: int _gr_poly_mulmod_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``.

    It is required that ``finv`` is the inverse of the reverse of
    ``f`` mod ``x^lenf``.

    Aliasing of ``res`` with any of the inputs is not permitted.

.. function:: int gr_poly_mulmod_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)

    Sets ``res`` to the remainder of the product of ``poly1``
    and ``poly2`` upon polynomial division by ``f``. ``finv``
    is the inverse of the reverse of ``f``.

.. function:: int _gr_poly_powmod_fmpz_binexp(gr_ptr res, gr_srcptr poly, const fmpz_t e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
              int gr_poly_powmod_fmpz_binexp(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, const gr_poly_t f, gr_ctx_t ctx)
              int _gr_poly_powmod_fmpz_binexp_preinv(gr_ptr res, gr_srcptr poly, const fmpz_t e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
              int gr_poly_powmod_fmpz_binexp_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
              int _gr_poly_powmod_x_fmpz_preinv(gr_ptr res, const fmpz_t e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
              int gr_poly_powmod_x_fmpz_preinv(gr_poly_t res, const fmpz_t e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
              int _gr_poly_powmod_ui_binexp(gr_ptr res, gr_srcptr poly, ulong e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
              int gr_poly_powmod_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong e, const gr_poly_t f, gr_ctx_t ctx)
              int _gr_poly_powmod_ui_binexp_preinv(gr_ptr res, gr_srcptr poly, ulong e, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
              int gr_poly_powmod_ui_binexp_preinv(gr_poly_t res, const gr_poly_t poly, ulong e, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
              int _gr_poly_powmod_fmpz_sliding_preinv(gr_ptr res, gr_srcptr poly, const fmpz_t e, ulong k, gr_srcptr f, slong lenf, gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
              int gr_poly_powmod_fmpz_sliding_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e, ulong k, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)

    Sets ``res`` to ``poly`` raised to the power ``e``
    modulo ``f``.

.. function:: int _gr_poly_compose_mod_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
              int gr_poly_compose_mod_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
              int _gr_poly_compose_mod_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
              int gr_poly_compose_mod_brent_kung(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)
              int _gr_poly_compose_mod(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_ctx_t ctx)
              int gr_poly_compose_mod(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, gr_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h` where *f*, *g*
    and *h* are given by *poly1*, *poly2* and *poly3*.

    The underscore methods require that the lengths are nonzero
    and that the length of `g` is one less than the length of `h` (possibly
    with zero padding). We also require that the length of `f` is less than
    the length of `h`. The underscore methods do not support aliasing.

.. function:: int _gr_poly_compose_mod_horner_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
              int gr_poly_compose_mod_horner_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)
              int _gr_poly_compose_mod_brent_kung_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
              int gr_poly_compose_mod_brent_kung_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)
              int _gr_poly_compose_mod_preinv(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, gr_srcptr poly3, slong len3, gr_srcptr poly3inv, slong inv3len, gr_ctx_t ctx)
              int gr_poly_compose_mod_preinv(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, const gr_poly_t poly3, const gr_poly_t poly3inv, gr_ctx_t ctx)

    Versions accepting an additional precomputed argument ``poly3inv``
    representing the inverse of the reverse of ``poly3``.

.. function:: int _gr_poly_reduce_matrix_mod_poly(gr_mat_t A, const gr_mat_t B, const gr_poly_t f, gr_ctx_t ctx)

    Sets the ith row of ``A`` to the reduction of the ith row of `B` modulo
    `f` for `i=1,\ldots,\sqrt{\deg(f)}`. We require `B` to be at least
    a `\sqrt{\deg(f)}\times \deg(f)` matrix and `f` to be nonzero.

.. function:: int _gr_poly_precompute_matrix(gr_mat_t A, gr_srcptr f, gr_srcptr g, slong leng, gr_srcptr ginv, slong lenginv, gr_ctx_t ctx)

    Sets the ith row of ``A`` to `f^i` modulo `g` for
    `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    be the inverse of the reverse of ``g`` and `g` to be nonzero.

.. function:: int gr_poly_precompute_matrix(gr_mat_t A, const gr_poly_t f, const gr_poly_t g, const gr_poly_t ginv, gr_ctx_t ctx)

    Sets the ith row of ``A`` to `f^i` modulo `g` for
    `i=1,\ldots,\sqrt{\deg(g)}`. We require `A` to be a
    `\sqrt{\deg(g)}\times \deg(g)` matrix. We require ``ginv`` to
    be the inverse of the reverse of ``g``.

.. function:: int _gr_poly_compose_mod_brent_kung_precomp_preinv(gr_ptr res, gr_srcptr f, slong lenf, const gr_mat_t A, gr_srcptr h, slong lenh, gr_srcptr hinv, slong lenhinv, gr_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that `h` is nonzero. We require that the ith row of `A` contains
    `g^i` for `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a
    `\sqrt{\deg(h)}\times \deg(h)` matrix. We also require that the
    length of `f` is less than the length of `h`. Furthermore, we
    require ``hinv`` to be the inverse of the reverse of ``h``.
    The output is not allowed to be aliased with any of the inputs.

    The algorithm used is the Brent-Kung matrix algorithm.

.. function:: int gr_poly_compose_mod_brent_kung_precomp_preinv(gr_poly_t res, const gr_poly_t f, const gr_mat_t A, const gr_poly_t h, const gr_poly_t hinv, gr_ctx_t ctx)

    Sets ``res`` to the composition `f(g)` modulo `h`. We require
    that the ith row of `A` contains `g^i` for
    `i=1,\ldots,\sqrt{\deg(h)}`, i.e. `A` is a `\sqrt{\deg(h)}\times
    \deg(h)` matrix. We require that `h` is nonzero and that `f` has
    smaller degree than `h`. Furthermore, we require ``hinv`` to be
    the inverse of the reverse of ``h``. This version of Brent-Kung
    modular composition is particularly useful if one has to perform
    several modular composition of the form `f(g)` modulo `h` for
    fixed `g` and `h`.


Test functions
-------------------------------------------------------------------------------

The following functions run *iters* test iterations, generating
polynomials up to length *maxn*. If *ctx* is set to ``NULL``, a random
ring is generated on each test iteration, otherwise the given ring is used.

.. function:: void _gr_poly_test_mullow(gr_method_poly_binary_trunc_op mullow_impl, gr_method_poly_binary_trunc_op mullow_ref, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``mullow_impl`` for correctness as an implementation
    of :func:`_gr_poly_mullow`. 
    A reference implementation to compare against can be provided as
    ``mullow_ref``; if ``NULL``, classical multiplication is used.

.. function:: void _gr_poly_test_divrem(gr_method_poly_binary_binary_op divrem_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``divrem_impl`` for correctness as an implementation
    of :func:`_gr_poly_divrem`.

.. function:: void _gr_poly_test_div(gr_method_poly_binary_op div_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``div_impl`` for correctness as an implementation
    of :func:`_gr_poly_div`.

.. function:: void _gr_poly_test_inv_series(gr_method_poly_unary_trunc_op inv_series_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``inv_series_impl`` for correctness as an implementation
    of :func:`_gr_poly_inv_series`.

.. function:: void _gr_poly_test_div_series(gr_method_poly_binary_trunc_op div_series_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``div_series_impl`` for correctness as an implementation
    of :func:`_gr_poly_div_series`.

.. function:: void _gr_poly_test_gcd(gr_method_poly_gcd_op gcd_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``gcd_impl`` for correctness as an implementation
    of :func:`_gr_poly_gcd`.

.. function:: void _gr_poly_test_xgcd(gr_method_poly_xgcd_op xgcd_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

    Tests the given function ``xgcd_impl`` for correctness as an implementation
    of :func:`_gr_poly_xgcd`.


.. raw:: latex

    \newpage

