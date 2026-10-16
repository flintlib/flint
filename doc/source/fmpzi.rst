.. _fmpzi:

**fmpzi.h** -- Gaussian integers
===============================================================================

This module allows working with elements of the ring `\mathbb{Z}[i]`.
At present, only a minimal interface is provided.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpzi_struct

.. type:: fmpzi_t

    Contains a pairs of integers representing the real and imaginary
    parts. An *fmpzi_t* is defined as an array
    of length one of type *fmpzi_struct*, permitting an *fmpzi_t* to
    be passed by reference.

.. macro:: fmpzi_realref(x)

    Macro giving a pointer to the real part of *x*.

.. macro:: fmpzi_imagref(x)

    Macro giving a pointer to the imaginary part of *x*.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void fmpzi_init(fmpzi_t x)

.. function:: void fmpzi_clear(fmpzi_t x)

.. function:: void fmpzi_swap(fmpzi_t x, fmpzi_t y)

.. function:: void fmpzi_zero(fmpzi_t x)

.. function:: void fmpzi_one(fmpzi_t x)

.. function:: void fmpzi_set(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_set_si_si(fmpzi_t res, slong a, slong b)

Input and output
-------------------------------------------------------------------------------

.. function:: void fmpzi_print(const fmpzi_t x)

Random number generation
-------------------------------------------------------------------------------

.. function:: void fmpzi_randtest(fmpzi_t res, flint_rand_t state, flint_bitcnt_t bits)

Properties
-------------------------------------------------------------------------------

.. function:: int fmpzi_equal(const fmpzi_t x, const fmpzi_t y)

.. function:: int fmpzi_is_zero(const fmpzi_t x)

.. function:: int fmpzi_is_one(const fmpzi_t x)

Units
-------------------------------------------------------------------------------

.. function:: int fmpzi_is_unit(const fmpzi_t x)

.. function:: slong fmpzi_canonical_unit_i_pow(const fmpzi_t x)

.. function:: void fmpzi_canonicalise_unit(fmpzi_t res, const fmpzi_t x)

Norms
-------------------------------------------------------------------------------

.. function:: slong fmpzi_bits(const fmpzi_t x)

.. function:: void fmpzi_norm(fmpz_t res, const fmpzi_t x)

Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpzi_conj(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_neg(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sqr(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)

Division
-------------------------------------------------------------------------------

.. function:: void fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)

    Sets *q* to the quotient of *x* and *y*, assuming that the
    division is exact.

.. function:: void fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

    Computes a quotient and remainder satisfying
    `x = q y + r` with `N(r) \le N(y)/2`, with a canonical
    choice of remainder when breaking ties.

.. function:: void fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)

    Computes a quotient and remainder satisfying
    `x = q y + r` with `N(r) < N(y)`, with an implementation-defined,
    non-canonical choice of remainder.

.. function:: slong fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x)

    Divide *x* exactly by the largest possible power `(1+i)^k`
    and return the exponent *k*.

GCD
-------------------------------------------------------------------------------

.. function:: void fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd_shortest(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
              void fmpzi_gcd(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

    Computes the GCD of *x* and *y*. The result is in canonical
    unit form.

    The *euclidean* version is a straightforward implementation
    of Euclid's algorithm. The *euclidean_improved* version is
    optimized by performing approximate divisions.
    The *binary* version uses a (1+i)-ary analog of the binary
    GCD algorithm for integers [Wei2000]_.
    The *shortest* version finds the GCD as the shortest vector in a lattice.
    The default version chooses an algorithm automatically.

Primality testing
--------------------------------------------------------------------------------

.. function:: int fmpzi_is_prime(const fmpzi_t n)

    Check whether `n` is a Gaussian prime.

.. function:: int fmpzi_is_probabprime(const fmpzi_t n)

    Check whether `n` is a probable Gaussian prime.

.. function:: int _fmpzi_poly_mulmid_classical_fft_small(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

    Computes the coefficients `[\mathtt{nlo}, \mathtt{nhi})` of the
    product of the Gaussian integer polynomials ``(poly1, len1)`` and
    ``(poly2, len2)``, with the coefficient arithmetic done pointwise
    on fft_small transforms: two transforms per coefficient (real and
    imaginary parts) and four pointwise operations per complex
    product, the shorter operand's transforms kept and the longer
    rolling as in :func:`_fmpz_poly_mulmid_classical_fft_small`.
    Squaring shares transforms, pairs symmetric terms and computes
    diagonal squares with three multiplications. The lengths may be
    given in either order. Returns 1 on success and 0 when the
    transformed representation is unavailable or unprofitable in which
    case nothing is written and the caller should fall back to a
    different routine.

.. function:: void _fmpzi_poly_mulmid(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

    Sets ``(res, nhi - nlo)`` to the coefficients at indices
    `[\mathtt{nlo}, \mathtt{nhi})` of the product of the Gaussian integer
    polynomials ``(poly1, len1)`` and ``(poly2, len2)``, choosing an
    algorithm. Assumes that ``len1`` and ``len2`` are positive and that
    `0 \le nlo < nhi \le len1 + len2 - 1`; the lengths may be given in
    either order, and ``res`` must not alias the inputs. This is the
    multiplication used by the ``fmpzi`` generic ring for polynomial
    multiplication, low products and middle products.

    The choice depends on `b`, the average over the two operands of the
    larger part size in bits, and on `q`, the number of coefficient
    products in the window per operand coefficient (`n/2` for an
    `n \times n` product, `n/3` for the Newton middle product
    `(2n - 1) \times n \to [n - 1, 2n - 1)`): the classical algorithm
    for windows of at most 12 products; four real products when a part
    vanishes or is much smaller than the other (they skip the vanishing
    products); :func:`_fmpzi_poly_mulmid_classical_fft_small` from 8192
    bits up to a bound on `q` increasing with `b`; real products, three
    (Karatsuba) from about 256 bits, below a bound on `q` decreasing with
    `b`; and otherwise :func:`_fmpzi_poly_mulmid_KS`, except for
    :func:`_fmpzi_poly_mulmid_SS` on squares at 384-3072 bits, short full
    products at 768-1536 bits and middle products whose cyclic transform
    saves a quarter of the length at 768-3072 bits (short ones) and from
    24576 bits, and for huge products (`m_1 + m_2 \ge 32768` and
    `(len_1 + len_2)(m_1 + m_2) > 2.5 \cdot 10^8` for the larger part sizes
    `m_1, m_2` of the operands), where Schoenhage-Strassen runs within
    about 10% of the time in about a fifth of the memory. The boundaries were
    fitted to timings on x86-64 with uniformly random coefficients.

.. function:: void _fmpzi_poly_mulmid_KS(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

    Sets ``(res, nhi - nlo)`` to the coefficients at indices
    `[\mathtt{nlo}, \mathtt{nhi})` of the product of the Gaussian integer
    polynomials ``(poly1, len1)`` and ``(poly2, len2)``, using Kronecker
    substitution. Assumes that ``len1`` and ``len2`` are positive and
    that `0 \le nlo < nhi \le len1 + len2 - 1`; the lengths may be given
    in either order, and ``res`` may alias the inputs.

    The real and imaginary parts of each operand are packed into
    integers at a common field size, and a single
    :func:`flint_mpn_mul_complex` (or :func:`flint_mpn_sqr_complex`)
    call forms both parts of the packed product. Above the transform
    cutoff of those functions the four packed parts are thus
    transformed once and combined pointwise (four forward and two
    inverse transforms, two and two for squaring), and below it their
    Karatsuba or schoolbook schemes apply. Windows are computed as
    two's-complement windows of the packed parts, which unpack directly
    into signed coefficients: with :func:`flint_mpn_mullow_complex` (or
    :func:`flint_mpn_sqrlow_complex`) when the window starts at the
    bottom and ends below the top of the product, and with
    :func:`flint_mpn_mulmid_complex` (or :func:`flint_mpn_sqrmid_complex`)
    when it starts above the bottom, three guard limbs below the first
    field absorbing the error of the windowed products. Above the
    transform cutoff the latter reconstructs only the CRT coefficients
    reaching the window and, when the window allows it, uses a cyclic
    transform, which is where a middle product saves transform length. When one operand is real,
    the products are formed separately, without sharing the transform
    of the other operand.

.. function:: void _fmpzi_poly_mulmid_SS(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)
              void _fmpzi_poly_mulmid_SS_pair(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)
              void _fmpzi_poly_mulmid_SS_single(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1, const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)

    Sets ``(res, nhi - nlo)`` to the coefficients at indices
    `[\mathtt{nlo}, \mathtt{nhi})` of the product of the Gaussian integer
    polynomials ``(poly1, len1)`` and ``(poly2, len2)``, using
    Schoenhage-Strassen convolution. Assumes that ``len1`` and ``len2``
    are positive and that `0 \le nlo < nhi \le len1 + len2 - 1`; the
    lengths may be given in either order, and ``res`` may alias the
    inputs.

    In the coefficient ring `\mathbb{Z}/(2^N+1)` of the transform,
    `s = 2^{N/2}` is a square root of `-1`, and
    `\mathbb{Z}[i]/(2^N+1) \cong \mathbb{Z}/(2^N+1) \times \mathbb{Z}/(2^N+1)`
    via `a + bi \mapsto (a + sb, a - sb)`. The *pair* method
    (``_SS_pair``) therefore computes the Gaussian product as two
    ordinary convolutions, of
    `u = \operatorname{re} + s \operatorname{im}` and
    `v = \operatorname{re} - s \operatorname{im}`, with a ring size
    covering the components of the result, and recovered as
    `\operatorname{re} = (U + V)/2`,
    `\operatorname{im} = (V - U) s / 2`; all conversions are shifts and
    additions. This uses six transforms and two pointwise products per
    point (four transforms for squaring). The *single embedding*
    (``_SS_single``) uses the homomorphism `a + bi \mapsto a + sb` alone:
    one convolution in a ring of twice the size gives
    `T = \operatorname{re} + s \operatorname{im}`, from whose two's
    complement representation both parts are read off, the ring being
    large enough that `|T| < 2^{N-1}`. The transform work and the memory
    are about the same as for the pair method, but the ring is rounded
    up to a multiple of the transform length once rather than twice,
    which favors it for small coefficients relative to the length; it
    also halves the number of pointwise products, which favors it for
    small rings, while for medium rings the pointwise products of the
    larger ring cost more than twice as much. The main function chooses
    between them by the rounded ring sizes, as measured on x86-64. The
    convolution
    length is
    chosen as by :func:`fft_convolution_choose_mid`, so middle products
    use the cyclic (wraparound) transform when it is shorter. When an
    operand is real or purely imaginary, its nonzero part is transformed
    once and multiplied by both parts of the other operand (two
    convolutions sharing that transform, five transforms in all), or by
    the single nonzero part of the other operand (one convolution).
