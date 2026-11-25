/*
    Copyright (C) 2025 Albin Ahlbäck
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_IMPL_H
#define NMOD_POLY_IMPL_H

#include "nmod_types.h"
#include "ulong_extras.h"  /* n_mulhi */
#include "nmod_vec.h"
#include "flint-mparam.h"  /* FLINT_NMOD_POLY_IFMA_* */
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"    /* _nmod_poly_should_directly_fft */
#endif

/*
    Products with a short operand by rows (_nmod_poly_mulmid_rows): one
    _nmod_vec_scalar_addmul_nmod of the longer operand per coefficient of
    the shorter one, of length len2 from 3 to _NMOD_POLY_ROWS_MAX_LEN2,
    rather than one dot product of at most len2 terms per output coefficient
    as in the classical code. _nmod_poly_mul, _nmod_poly_mullow and
    _nmod_poly_mulmid use them in place of the classical code and of KS when
    these scalar products have a vector kernel for the modulus and every row
    is long enough: the shortest row, of length row_len, reaches the length
    of the kernels and ratio * len2, since a row costs a few times more than
    a dot product on top of their common arithmetic. For powers of 2, whose
    dot products need no reduction, only up to 2^31 (2^48 with the 512-bit
    kernels) and with rows twice as long, 2 ratio len2: there the rows are
    1.1-4x faster than the classical code, which wins on shorter rows and
    from 2^32 with AVX2 (the fp50 kernel) or 2^50 with AVX-512 (measured Nov
    2026).

    _NMOD_POLY_ROWS_RATIO: the ratio for full and low products, 4 with the
    AVX2 kernels above 32 bits (doubles, peeling to an aligned address,
    scalar tails), 2 with the 512-bit kernels (the only ones above 50 bits,
    masked tails) and up to 32 bits.

    _NMOD_POLY_ROWS_MAX_LEN2: the largest len2, 12 with the NEON kernels
    (on the M4, at 31-32 bits, where the classical code does not
    accumulate in one word, the rows are 1.2-1.6x faster up to len2 =
    12-13, measured Nov 2026).

    _NMOD_POLY_NEON: 1 with the NEON kernels (aarch64), where KS compares
    better with the classical code in one word (on the M4, see mul.c,
    mullow.c and _nmod_poly_mulmid_want_one_word).
*/
#define _NMOD_POLY_ROWS_RATIO(mod) \
    ((NMOD_VEC_SIMD_SCALAR_MAX_BITS(1) > 50 || NMOD_BITS(mod) <= 32) ? 2 : 4)

#if defined(__ARM_NEON) && defined(__aarch64__)
# define _NMOD_POLY_NEON 1
# define _NMOD_POLY_ROWS_MAX_LEN2 12
#else
# define _NMOD_POLY_NEON 0
# define _NMOD_POLY_ROWS_MAX_LEN2 16
#endif

static inline int
_nmod_poly_want_rows(slong row_len, slong len2, slong ratio, nmod_t mod)
{
    return len2 >= 3 && len2 <= _NMOD_POLY_ROWS_MAX_LEN2
        && NMOD_BITS(mod) <= NMOD_VEC_SIMD_SCALAR_MAX_BITS(1)
        && row_len >= FLINT_MAX(NMOD_VEC_SIMD_MIN_LEN, ratio * len2)
        && ((mod.n & (mod.n - 1)) != 0
            || (NMOD_BITS(mod) <= (NMOD_VEC_SIMD_SCALAR_MAX_BITS(1) > 50 ? 49 : 32)
                && row_len >= 2 * ratio * len2));
}

void _nmod_poly_mulmid_rows(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

/* whether sums of len2 products of coefficients fit in one word: (n - 1)^2
   len2 < 2^64, which implies n <= 2^32 for len2 >= 2 */
static inline int
_nmod_poly_one_word(slong len2, nmod_t mod)
{
    return NMOD_BITS(mod) <= FLINT_BITS / 2
        && n_mulhi((mod.n - 1) * (mod.n - 1), len2) == 0;
}

/*
    Without NEON, the classical code accumulating in one word beats KS (or
    KS2, KS4) for a shorter operand of len2 < c(b) coefficients, b the
    number of bits of n, when the longer one has at least 2 len2 (full
    products, _nmod_poly_mul) or is cut to n >= 2 len2 (low products): c(b)
    = max(12, 2b) below 22 bits, 4b - 24 from 22 bits (fft_small takes
    over at cutoff_len 80-100 up to 21 bits, about 300 above). Measured
    crossovers on a Cascade Lake, gcc / clang, len1 >= 4 len2 (Nov 2026):
    12 / 9 at 4 bits, 16 / 12 at 8 bits, 18-32 / 18-20 at 12 bits, up to
    fft_small at 18-20 bits with gcc, 72-112 / 36-96 at 22-28 bits; Ice
    Lake and Zen 4 (without their IFMA kernel) and low products n x n/4 at
    least as far.
*/
static inline slong
_nmod_poly_one_word_cutoff(slong bits)
{
    return (bits < 22) ? FLINT_MAX(12, 2 * bits) : 4 * bits - 24;
}

/*
    Middle products by rows accumulating in one word
    (_nmod_poly_mulmid_rows_one_word): _nmod_poly_mulmid uses them for m
    output coefficients and len2 = 3 to 16 below 6 bits, 32 below 10 bits
    and 64 above, when m >= max(8, len2), in place of the rows (from len2 =
    6; below, only where the rows do not apply), of the classical code (one
    dot product per output, which has no one-word path), of KS and of
    fft_small; not for full products, which go to _nmod_poly_mul, nor for
    squares, where the classical code computes half of the products. With
    NEON, for m >= 16 and len2 up to 10 below 14 bits, 64 above, and before
    the rows. Measured (Nov 2026): within 3% of the best of the other
    algorithms from 8 outputs for len2 = 5 to 14 with gcc, 10 to 16 with
    clang, from about len2 outputs above (Cascade Lake, Ice Lake, Zen 4);
    on the M4 from 16 outputs, and KS faster for len2 >= 12 below 14 bits.
*/
static inline int
_nmod_poly_mulmid_want_one_word(slong len2, slong m, nmod_t mod)
{
    slong bits = NMOD_BITS(mod);

    if (_NMOD_POLY_NEON)
        return m >= 16 && len2 >= 3 && len2 <= (bits < 14 ? 10 : 64)
            && _nmod_poly_one_word(len2, mod);

    return m >= FLINT_MAX(8, len2) && len2 >= 3
        && len2 <= (bits < 6 ? 16 : bits < 10 ? 32 : 64)
        && _nmod_poly_one_word(len2, mod);
}

void _nmod_poly_mulmid_rows_one_word(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

/* the shortest row of _nmod_poly_mulmid_rows: the row length is a concave
   function of the index of the coefficient of poly2, so it is that of the
   first or the last coefficient used */
static inline slong
_nmod_poly_mulmid_rows_min_len(slong len1, slong len2, slong nlo, slong nhi)
{
    slong j0 = FLINT_MAX(0, nlo - len1 + 1);
    slong j1 = FLINT_MIN(len2, nhi) - 1;

    return FLINT_MIN(FLINT_MIN(j0 + len1, nhi) - FLINT_MAX(j0, nlo),
                     FLINT_MIN(j1 + len1, nhi) - FLINT_MAX(j1, nlo));
}

/*
    Classical products with AVX512-IFMA (_nmod_poly_mulmid_ifma, any window
    of the product, any n): for a shorter operand of
    _NMOD_POLY_IFMA_MIN_LEN(bits) to at most 511 coefficients, used in place
    of the rows, of the classical code, of KS and of fft_small. The
    function exists in every build (it runs the classical code without
    IFMA), the dispatchers only call it when NMOD_POLY_HAVE_IFMA.
    _nmod_poly_mulmid_ifma_crt is the same, but takes the variant by CRT
    for n > 2^52 whatever the lengths (for the tests), and
    _nmod_poly_mulmid_ifma_ptr takes mod by pointer (for the dispatchers,
    see below).

    _NMOD_POLY_IFMA_MIN_LEN: below 8 coefficients the rows are faster
    (rows / kernel 0.82-0.86 at 7 on Ice Lake, 31 to 50 bits), except for
    52-bit moduli, where the rows run the slower scalar products of 51 to
    63 bits: the kernel from 4 (1.05-1.22x faster at 4, 1.5-1.9x at 6 and
    7; Zen 4 and Ice Lake, Nov 2026). Above 52 bits the kernel has larger
    fixed costs (the shorter operand split or reduced, the two-limb REDC):
    it starts from FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52 coefficients for
    53- to 63-bit moduli (Ice Lake 8, at least 0.99x the dispatchers
    without it; Zen 4 12: 0.93-0.95x at 10, 0.99-1.00x at 12), from 4 for
    64-bit moduli, which have no rows (1.5x or more at 4 on both), and is
    used for odd n and powers of 2 (for even n, the scalar reduction of
    each output, NMOD_RED3, makes it 0.1-0.75x the dispatchers without it
    up to 16 coefficients), from 16 output coefficients (nout) and for the
    shorter shapes of _nmod_poly_ifma_above_52_in_range.

    The largest len2 (above, fft_small) comes from flint-mparam.h:
    FLINT_NMOD_POLY_IFMA_MAX_LEN(bits) for full and middle products, the
    crossover with the dispatchers without the kernel on products of
    lengths L and L, 4L and L, 16L and L (the unbalanced ones cross first).
    It follows the number of 50-bit primes fft_small needs: one up to 20
    bits, two up to 45 bits (the kernel takes one IFMA per product up to 26
    bits, two above), three above (seven IFMA above 52 bits, or the CRT);
    Zen 4 crosses about twice as early as Ice Lake from 21 bits on. Modulo
    the 50-bit primes of fft_small, which need no CRT,
    FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME. Low products cross later: up to
    511 coefficients (at the primes of fft_small,
    FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME). Above 52 bits, the CRT
    variant from FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN (0: never): on Zen 4,
    where vpmullq is fast, it is 1.25-1.5x faster than seven IFMA from 64
    coefficients; on Ice Lake, 1.3-2x slower. The bound above 52 bits is
    256 on Ice Lake (448 for balanced products) and 128 on Zen 4 (with the
    CRT; 192 for balanced products).

    For balanced shapes (full products of len1 <= 2 len2 + 1, middle
    products of at most 3 len2 outputs), the bound is
    FLINT_NMOD_POLY_IFMA_BALANCED_MAX_LEN(bits), at least the one above (on
    Zen 4, 256 instead of 192 at 27-45 bits), but beyond
    FLINT_NMOD_POLY_IFMA_MAX_LEN(bits) not modulo the primes that fft_small
    transforms directly (_nmod_poly_fft_small_direct).

    With the one-word Karatsuba (see below), from flint-mparam.h: full and
    middle products of blocks (the shorter operand, or the outputs of a
    middle product if fewer) longer than FLINT_NMOD_POLY_IFMA_KARA_MAX_LEN
    go to the Karatsuba, which is faster there. Low products only from
    blocks of FLINT_NMOD_POLY_IFMA_KARA_LOW_MIN_LEN, unbalanced, or balanced
    from 4 times that: the kernel computes only the terms of the outputs.
    At 27-30 bits, where the kernel takes two IFMA per product, one block of
    FLINT_NMOD_POLY_IFMA_KARA_BLOCK_MIN_LEN to 16 coefficients too, for
    balanced full products and squares, and middle products by 16
    coefficients of fewer than 128 outputs (if that bound is at most 16).
*/
#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512DQ__) \
        && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_POLY_HAVE_IFMA 1
#else
# define NMOD_POLY_HAVE_IFMA 0
#endif

#if !defined(FLINT_NMOD_POLY_IFMA_MAX_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME) \
    || !defined(FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME) \
    || !defined(FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52) \
    || !defined(FLINT_NMOD_POLY_IFMA_BALANCED_MAX_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_KARA_MAX_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_KARA_BLOCK_MIN_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_KARA_LOW_MIN_LEN) \
    || !defined(FLINT_NMOD_POLY_KARA_AVX512_RULES) \
    || !defined(FLINT_NMOD_POLY_KARA_ONE_PRIME_MAX_LEN) \
    || !defined(FLINT_NMOD_POLY_KARA_TWO_PRIMES_MAX_LEN)
# error "flint-mparam.h must define FLINT_NMOD_POLY_IFMA_MAX_LEN, FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME, FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME, FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN, FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52, FLINT_NMOD_POLY_IFMA_BALANCED_MAX_LEN, FLINT_NMOD_POLY_IFMA_KARA_MAX_LEN, FLINT_NMOD_POLY_IFMA_KARA_BLOCK_MIN_LEN, FLINT_NMOD_POLY_IFMA_KARA_LOW_MIN_LEN, FLINT_NMOD_POLY_KARA_AVX512_RULES, FLINT_NMOD_POLY_KARA_ONE_PRIME_MAX_LEN and FLINT_NMOD_POLY_KARA_TWO_PRIMES_MAX_LEN"
#endif

#define _NMOD_POLY_IFMA_MIN_LEN(bits) \
    ((bits) < 52 ? 8 : (bits) == 52 || (bits) == 64 ? 4 \
     : FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52)

/* whether n is one of the primes of the default fft_small context */
int _nmod_poly_is_fft_small_prime(ulong n);

/*
    Whether fft_small computes the product of operands of lengths
    len1 >= len2 modulo the prime n itself, with one transform and no CRT
    (fft_small's own test, with its cache of the last prime and the last
    composite; the depth is that of a full product).
*/
static inline int
_nmod_poly_fft_small_direct(slong len1, slong len2, ulong n)
{
#if FLINT_HAVE_FFT_SMALL
    return _nmod_poly_should_directly_fft(len1,
                                    FLINT_BIT_COUNT(len1 + len2 - 2), n);
#else
    (void) len1;
    (void) len2;
    (void) n;
    return 0;
#endif
}

/*
    Above 52 bits, the shapes where the kernel takes odd n and powers of 2:
    from 16 output coefficients; middle products from 8 outputs of a
    shorter operand of at least 12; low products from 14; from 62 bits,
    every kind from 11 (low products from 12). They depend on the lengths
    and the bit size only: the dispatchers send the other shapes to the
    other algorithms at once.
*/
static inline int
_nmod_poly_ifma_above_52_in_range(slong len2, slong nout, slong bits,
                                  int low, int mid)
{
    return nout >= 16 || (mid && nout >= 8 && len2 >= 12)
        || (low && nout >= 14) || (bits >= 62 && nout >= (low ? 12 : 11));
}

/* nout: the number of output coefficients; low: whether the product is a
   low product (mullow). _nmod_poly_ifma_in_range: the conditions on the
   lengths, which depend on the bit size of n only (see the dispatchers).
   Low products and middle products of few outputs (len2 >= 4 nout) go up
   to the kernel's 511: for the latter it costs nout len2 products, and for
   48-96 outputs by 4 nout + r it was 1.25-5.4x faster than the Karatsuba
   or fft_small on Zen 4 and Ice Lake (Nov 2026) */
static inline int
_nmod_poly_ifma_in_range(slong len2, slong nout, slong bits, int low)
{
#if NMOD_POLY_HAVE_IFMA
    return len2 >= _NMOD_POLY_IFMA_MIN_LEN(bits)
        && len2 <= ((low || 4 * nout <= len2) ? 511 : (nout <= 3 * len2)
                    ? FLINT_NMOD_POLY_IFMA_BALANCED_MAX_LEN(bits)
                    : FLINT_NMOD_POLY_IFMA_MAX_LEN(bits));
#else
    (void) len2;
    (void) nout;
    (void) bits;
    (void) low;
    return 0;
#endif
}

static inline int
_nmod_poly_want_ifma(slong len1, slong len2, slong nout, nmod_t mod, int low)
{
#if NMOD_POLY_HAVE_IFMA
    slong bits = NMOD_BITS(mod);

    if (!_nmod_poly_ifma_in_range(len2, nout, bits, low))
        return 0;

    /* the longer balanced shapes (but middle products of fewer outputs
       than len2): not where fft_small transforms n directly */
    if (!low && nout >= len2 && len2 > FLINT_NMOD_POLY_IFMA_MAX_LEN(bits)
            && _nmod_poly_fft_small_direct(len1, len2, mod.n))
        return 0;

    /* bits > 52 rather than n > 2^52, as in _nmod_poly_dispatch: for
       n = 2^52 the kernel takes its 52-bit path anyway */
    if (bits > 52)
        return _nmod_poly_ifma_above_52_in_range(len2, nout, bits, low,
                                        !low && nout < len1 + len2 - 1)
            && ((mod.n & 1) || (mod.n & (mod.n - 1)) == 0);

    return len2 <= (low ? FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME
                        : FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME)
        || bits != 50 || !_nmod_poly_is_fft_small_prime(mod.n);
#else
    /* no kernel (and 2^52 does not fit in a 32-bit word) */
    (void) len1;
    (void) len2;
    (void) nout;
    (void) mod;
    (void) low;
    return 0;
#endif
}

void _nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);
void _nmod_poly_mulmid_ifma_crt(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);
void _nmod_poly_mulmid_ifma_ptr(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, const nmod_t * mod);

/*
    Karatsuba with AVX2 (_nmod_poly_mulmid_karatsuba, any window of the
    product, see mulmid_karatsuba.c). Two arithmetics: one word where its
    sums fit (_nmod_poly_karatsuba_one_word), else doubles for n < 2^50
    (with FMA). The function exists in every build: without AVX2, or for
    n >= 2^50, it runs the classical code. _nmod_poly_mulmid_karatsuba_fp50
    forces the doubles (for the tests).
*/
#if FLINT_BITS == 64 && defined(__AVX2__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_POLY_HAVE_KARATSUBA 1
#else
# define NMOD_POLY_HAVE_KARATSUBA 0
#endif

/* MSVC does not define __FMA__ but /arch:AVX2 has FMA (as in
   machine_vectors.h) */
#if NMOD_POLY_HAVE_KARATSUBA && (defined(__FMA__) || defined(_MSC_VER))
# define NMOD_POLY_HAVE_KARATSUBA_FP50 1
#else
# define NMOD_POLY_HAVE_KARATSUBA_FP50 0
#endif

void _nmod_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);
void _nmod_poly_mulmid_karatsuba_fp50(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

/*
    Whether the one-word arithmetic is exact for blocks of len coefficients
    (padded to a multiple of 16) and output coefficients of at most terms
    products. After d halvings of the blocks down to 16 coefficients, the
    operands are sums of up to 2^d coefficients, which must stay below 2^32
    (2^31 in absolute value for the signed middle products); the output
    coefficients, at most terms (n - 1)^2, must stay below 2^64. This gives
    n up to 2^26 for 1024 coefficients, 2^28 for 256 and 2^30 for 16 (one
    bit less for middle products). _nmod_poly_karatsuba_one_word_bits:
    whether this holds for every n of the given bit size.
*/
static inline int
_nmod_poly_karatsuba_one_word(ulong n, slong len, slong terms, int mid)
{
    /* the number of halvings of the (len + 15) / 16 blocks */
    slong d = FLINT_BIT_COUNT((len + 15) / 16 - 1);

    return n <= ((UWORD(1) << (mid ? 31 : 32)) >> d)
        && n_mulhi((n - 1) * (n - 1), (ulong) terms) == 0;
}

static inline int
_nmod_poly_karatsuba_one_word_bits(slong bits, slong len, slong terms, int mid)
{
    slong d = FLINT_BIT_COUNT((len + 15) / 16 - 1);

    return bits <= (mid ? 31 : 32) - d && terms <= (WORD(1) << (64 - 2 * bits));
}

/*
    The kinds of products: full products, squares, low products (nlo = 0)
    and middle products. The rules below describe a product of operands of
    lengths len1 >= len2 (after the trimming of a middle product) with m
    output coefficients by its block length s (len2, or min(len2, m) for a
    middle product) and the other dimension l (len1, m for a low product,
    max(len2, m) for a middle product); the shape is unbalanced if l >= 2s.
*/
#define _NMOD_POLY_KIND_MUL 0
#define _NMOD_POLY_KIND_SQR 1
#define _NMOD_POLY_KIND_LOW 2
#define _NMOD_POLY_KIND_MID 3

/* the kind of the product of coefficients nlo to nhi (squares: the same
   operand twice, as in _nmod_poly_mul) */
FLINT_FORCE_INLINE int
_nmod_poly_product_kind(nn_srcptr poly1, slong len1, nn_srcptr poly2,
                        slong len2, slong nlo, slong nhi)
{
    return (nlo > 0) ? _NMOD_POLY_KIND_MID
        : (nhi < len1 + len2 - 1) ? _NMOD_POLY_KIND_LOW
        : (poly1 == poly2 && len1 == len2) ? _NMOD_POLY_KIND_SQR
        : _NMOD_POLY_KIND_MUL;
}

/* the window [nlo, nhi) of a product: drop the coefficients of the
   operands that no output coefficient uses, and put the longer operand
   first (for _nmod_poly_mulmid and the Karatsuba) */
FLINT_FORCE_INLINE void
_nmod_poly_mulmid_trim(nn_srcptr * poly1, slong * len1, nn_srcptr * poly2,
                       slong * len2, slong * nlo, slong * nhi)
{
    slong t;

    *len1 = FLINT_MIN(*len1, *nhi);
    *len2 = FLINT_MIN(*len2, *nhi);

    t = *nlo - (*len2 - 1);
    if (t > 0)
    {
        *poly1 += t;
        *len1 -= t;
        *nlo -= t;
        *nhi -= t;
    }

    t = *nlo - (*len1 - 1);
    if (t > 0)
    {
        *poly2 += t;
        *len2 -= t;
        *nlo -= t;
        *nhi -= t;
    }

    if (*len1 < *len2)
    {
        FLINT_SWAP(nn_srcptr, *poly1, *poly2);
        FLINT_SWAP(slong, *len1, *len2);
    }
}

/* the rule of the Karatsuba (mulmid_karatsuba.c), up to 50 bits */
int _nmod_poly_want_karatsuba(slong len1, slong len2, slong m, int kind, ulong n);

/*
    The doubles from 8 coefficients rather than 12, for 33 to 50 bits,
    with the AVX2 rules and without IFMA: unbalanced full and low products
    and squares (at 8, the Karatsuba was 1.4-2.2x faster than the other
    algorithms on Comet Lake, Nov 2026). The dispatchers pass these shapes
    on at once.
*/
#if NMOD_POLY_HAVE_KARATSUBA_FP50 && !NMOD_POLY_HAVE_IFMA \
        && !FLINT_NMOD_POLY_KARA_AVX512_RULES
# define _NMOD_POLY_KARATSUBA_HAVE_FROM_8 1
#else
# define _NMOD_POLY_KARATSUBA_HAVE_FROM_8 0
#endif

static inline int
_nmod_poly_karatsuba_from_8(slong len2, slong m, int kind, slong bits)
{
    return _NMOD_POLY_KARATSUBA_HAVE_FROM_8 && len2 >= 8
        && bits > 32 && bits <= 50
        && ((kind == _NMOD_POLY_KIND_MUL && m >= 3 * len2 - 1)
            || (kind == _NMOD_POLY_KIND_LOW && m >= 2 * len2)
            || kind == _NMOD_POLY_KIND_SQR);
}

/*
    Necessary conditions for the rule, from the lengths and the bit size of
    n only: the dispatchers check them first, and the tests check that they
    never decline a product that the rule takes. They are the shortest
    lengths of the rules:
    - from 12 coefficients (8 for the shapes of _nmod_poly_karatsuba_from_8),
      up to 50 bits;
    - two blocks only from 20 coefficients for squares and low products,
      and above 30 bits for every kind;
    - squares from 13 below 12 bits, not at 20-23 above 32 bits;
    - at 31-32 bits, unbalanced shapes from 16 (24 with the AVX-512 rules),
      balanced low products not at 17-27;
    - above 32 bits, balanced low products at 16 and from 28 (only from 32
      with the AVX-512 rules), middle products of few outputs (len2 >= 4m)
      from 16.
*/
static inline int
_nmod_poly_karatsuba_in_range(slong len2, slong m, int kind, slong bits)
{
#if NMOD_POLY_HAVE_KARATSUBA
    slong s;
    int low = (kind == _NMOD_POLY_KIND_LOW), unbalanced;

    if (len2 < 12 || m < 12 || bits > 50)
        return _nmod_poly_karatsuba_from_8(len2, m, kind, bits);

    s = (kind == _NMOD_POLY_KIND_MID) ? FLINT_MIN(len2, m) : len2;

    if (s > 16 && s < 20 && (bits > 30 || kind == _NMOD_POLY_KIND_SQR
                             || low))
        return 0;

    if (kind == _NMOD_POLY_KIND_SQR)
        return (bits >= 12 || s >= 13) && (bits <= 32 || s < 20 || s >= 24);

    if (bits <= 30)
        return 1;

    unbalanced = (kind == _NMOD_POLY_KIND_MID) ? FLINT_MAX(len2, m) >= 2 * s
               : low ? m >= 2 * len2 : m >= 3 * len2 - 1;

    if (bits <= 32)
        return unbalanced ? s >= (FLINT_NMOD_POLY_KARA_AVX512_RULES ? 24 : 16)
                          : (!low || s < 17 || s >= 28);

    if (low && !unbalanced)
        return FLINT_NMOD_POLY_KARA_AVX512_RULES ? s >= 32
                                                 : (s == 16 || s >= 28);

    if (kind == _NMOD_POLY_KIND_MID && len2 >= 4 * m)
        return s >= 16;

    return 1;
#else
    (void) len2;
    (void) m;
    (void) kind;
    (void) bits;
    return 0;
#endif
}

/*
    With IFMA, whether the one-word Karatsuba comes before the kernel (see
    FLINT_NMOD_POLY_IFMA_KARA_* above), from the lengths and the bit size
    of n: _nmod_poly_karatsuba_one_block_bits for one block, which it takes
    at once (one word holds there up to 30 bits);
    _nmod_poly_karatsuba_before_ifma_bits for that block or longer blocks,
    where the rule of the one-word Karatsuba decides.
    _nmod_poly_karatsuba_before_ifma: the latter from n (and one word
    holds).
*/
static inline int
_nmod_poly_karatsuba_one_block_bits(slong len1, slong len2, slong m,
                                    int kind, slong bits)
{
#if NMOD_POLY_HAVE_IFMA && NMOD_POLY_HAVE_KARATSUBA
    slong s = (kind == _NMOD_POLY_KIND_MID) ? FLINT_MIN(len2, m) : len2;

    if (bits > 30 || s > 16)
        return 0;

    /* balanced full products and squares, middle products by 16
       coefficients of fewer than 128 outputs */
    return (kind == _NMOD_POLY_KIND_MID)
        ? (len2 == 16 && m < 8 * 16
           && FLINT_NMOD_POLY_IFMA_KARA_BLOCK_MIN_LEN(bits) <= 16)
        : (kind != _NMOD_POLY_KIND_LOW && len1 < 2 * len2
           && s >= FLINT_NMOD_POLY_IFMA_KARA_BLOCK_MIN_LEN(bits));
#else
    (void) len1;
    (void) len2;
    (void) m;
    (void) kind;
    (void) bits;
    return 0;
#endif
}

static inline int
_nmod_poly_karatsuba_before_ifma_bits(slong len1, slong len2, slong m,
                                      int kind, slong bits)
{
#if NMOD_POLY_HAVE_IFMA && NMOD_POLY_HAVE_KARATSUBA
    slong s = (kind == _NMOD_POLY_KIND_MID) ? FLINT_MIN(len2, m) : len2;
    slong lo;

    if (bits > 30)
        return 0;

    if (s <= 16)
        return _nmod_poly_karatsuba_one_block_bits(len1, len2, m, kind, bits);

    /* low products: unbalanced, or balanced from 4 times the bound */
    if (kind == _NMOD_POLY_KIND_LOW)
    {
        lo = FLINT_NMOD_POLY_IFMA_KARA_LOW_MIN_LEN(bits);
        return lo > 0 && s >= lo && (m >= 2 * s || s >= 4 * lo);
    }

    return s > FLINT_NMOD_POLY_IFMA_KARA_MAX_LEN(bits);
#else
    (void) len1;
    (void) len2;
    (void) m;
    (void) kind;
    (void) bits;
    return 0;
#endif
}

static inline int
_nmod_poly_karatsuba_before_ifma(slong len1, slong len2, slong m, int kind,
                                 nmod_t mod)
{
#if NMOD_POLY_HAVE_IFMA && NMOD_POLY_HAVE_KARATSUBA
    int mid = (kind == _NMOD_POLY_KIND_MID);

    return _nmod_poly_karatsuba_before_ifma_bits(len1, len2, m, kind,
                                                 NMOD_BITS(mod))
        && mod.n < (UWORD(1) << 30)
        && _nmod_poly_karatsuba_one_word(mod.n,
                                    mid ? FLINT_MIN(len2, m) : len2, len2, mid);
#else
    (void) len1;
    (void) len2;
    (void) m;
    (void) kind;
    (void) mod;
    return 0;
#endif
}

/*
    The padding check of the rules. Both engines pad s to b = ceil(s/16)
    blocks of 16 and split them in halves, so that 5, 7 or 9-11 blocks cost
    about as much as the next power of 2 or 3 times a power of 2, and
    nearly empty blocks waste their padding. A block count is filled if
    s >= 0.9 (16 b). One block passes, two from 20 coefficients; then
    - _LOOSE: anything but low products of 3 or 5 blocks not filled;
    - _STRICT: filled blocks, b in {3, 4, 6, 8, 12, 16, 24, 32} ({3, 4, 8}
      for squares);
    - _STRICT40: as _STRICT, but 3 blocks from 40 coefficients and 4 blocks
      pass, filled or not.
*/
#define _NMOD_POLY_KARATSUBA_PAD_LOOSE 0
#define _NMOD_POLY_KARATSUBA_PAD_STRICT 1
#define _NMOD_POLY_KARATSUBA_PAD_STRICT40 2

static inline int
_nmod_poly_karatsuba_blocks_ok(slong s, int kind, int pad)
{
    slong b = (s + 15) / 16;
    int full = (10 * s >= 9 * 16 * b);

    if (b <= 1)
        return 1;
    if (b == 2)
        return s >= 20;
    if (pad == _NMOD_POLY_KARATSUBA_PAD_LOOSE)
        return kind != _NMOD_POLY_KIND_LOW || full || (b != 3 && b != 5);
    if (pad == _NMOD_POLY_KARATSUBA_PAD_STRICT40 && (b == 3 || b == 4))
        return s >= 40;
    if (kind == _NMOD_POLY_KIND_SQR)
        return full && (b == 3 || b == 4 || b == 8);
    return full && (b == 3 || b == 4 || b == 6 || b == 8 || b == 12
                    || b == 16 || b == 24 || b == 32);
}

/*
    In one word, below 14 bits, where KS packs several coefficients per
    word and is fast, the Karatsuba takes mostly lengths that pad well:
    - squares below 12 bits only with filled blocks; at 12-13 bits not with
      a last block filled below 3/4 (from 33 coefficients), and with IFMA,
      beyond 192, only with the filled blocks of the products;
    - below 8 bits, beyond 64 coefficients only filled blocks (but middle
      products from 4 bits: any);
    - below 8 bits, squares, low products and unbalanced shapes only with
      one block, two filled blocks (28-32) or four (58-64, squares 60-64);
      from 6 bits also squares from 116, balanced low products of 97-128
      and unbalanced shapes with filled blocks.
    The rules for 2 to 7 bits are rough: there KS depends on the lengths
    and on the processor more than a rule on the bit size can follow.
*/
static inline int
_nmod_poly_karatsuba_small_moduli_ok(slong s, int kind, slong bits,
                                     int unbalanced)
{
    const int strict = _NMOD_POLY_KARATSUBA_PAD_STRICT;

    if (kind == _NMOD_POLY_KIND_SQR && bits >= 12)
        return (NMOD_POLY_HAVE_IFMA && s > 192)
            ? _nmod_poly_karatsuba_blocks_ok(s, _NMOD_POLY_KIND_MUL, strict)
            : !(s > 32 && 4 * s < 3 * 16 * ((s + 15) / 16));

    if (kind == _NMOD_POLY_KIND_SQR
            && !_nmod_poly_karatsuba_blocks_ok(s, kind, strict))
        return 0;

    if (bits >= 8)
        return 1;

    if (s > 64 && (kind != _NMOD_POLY_KIND_MID || bits < 4)
            && !_nmod_poly_karatsuba_blocks_ok(s, kind, strict))
        return 0;

    if (s <= 16 || kind == _NMOD_POLY_KIND_MID)
        return 1;

    if (kind == _NMOD_POLY_KIND_SQR)
        return (s >= 28 && s <= 32) || (s >= 60 && s <= 64)
            || (bits >= 6 && s >= 116);

    if (kind == _NMOD_POLY_KIND_LOW || unbalanced)
        return (s >= 28 && s <= 32) || (s >= 58 && s <= 64)
            || (bits >= 6 && (unbalanced
                    ? _nmod_poly_karatsuba_blocks_ok(s, kind, strict)
                    : (s >= 97 && s <= 128)));

    return 1;
}

/*
    In one word, the longest blocks when fft_small needs one prime, by kind
    and bit size: shorter below 14 bits, where KS is fast; from 12 bits,
    320 (fft_small steps up at 257). Without IFMA, when one prime holds by
    the sizes (2b + log2(len2) <= 50), FLINT_NMOD_POLY_KARA_ONE_PRIME_MAX_LEN
    from 12 bits: 320 by default; 544 on Broadwell and Skylake-family
    cores, past the next step of fft_small at 513 (but 448 for middle
    products of 449-512).
*/
static inline slong
_nmod_poly_karatsuba_one_prime_max(int kind, slong bits, slong s, int sizes)
{
    /* 0: the default top */
    static const short tab[7][4] = {
        /*                   mul   sqr   low   mid */
        /* bits < 4  */   {   64,   64,   64,  192 },
        /*      < 6  */   {  128,   64,   64,  192 },
        /*      < 8  */   {    0,  128,    0,  512 },
        /*      < 10 */   {    0,  128,    0,    0 },
        /*      < 12 */   {    0,  256,    0,    0 },
        /*      < 14 */   {    0,  256,    0,    0 },
        /*     >= 14 */   {    0,  192,    0,    0 }
    };
    slong r = (bits < 4) ? 0 : (bits < 14) ? (bits - 2) / 2 : 6;
    slong top = (bits < 12) ? 256 : 320;

#if !NMOD_POLY_HAVE_IFMA
    if (bits >= 12 && sizes)
    {
        top = FLINT_NMOD_POLY_KARA_ONE_PRIME_MAX_LEN;
        if (top > 448 && kind == _NMOD_POLY_KIND_MID && s > 448 && s <= 512)
            top = 448;
    }
#else
    (void) s;
    (void) sizes;
#endif

    return tab[r][kind] ? tab[r][kind] : top;
}

/*
    The rule in one word (n < 2^30, where _nmod_poly_karatsuba_one_word
    holds); the dispatchers call it with n = 0, and then it returns -1
    where the answer depends on n.
    - Shortest blocks: 12 (13 for squares below 12 bits); below 4 bits, only
      balanced full products and middle products, no squares.
    - Below 14 bits: _nmod_poly_karatsuba_small_moduli_ok.
    - Two blocks only from 20 for squares and low products.
    - With IFMA, up to the longest lengths of the kernel, only filled blocks
      (squares with the blocks of the products).
    - Longest blocks: where fft_small takes over. With one prime by the
      sizes, _nmod_poly_karatsuba_one_prime_max; for a prime n that it
      transforms directly (from 20 bits and 100 coefficients), 512 for
      every kind, since it sets that transform up at each product; with two
      primes, with IFMA 2048 for full products and 1024 for the others,
      without IFMA FLINT_NMOD_POLY_KARA_TWO_PRIMES_MAX_LEN (half of it for
      squares). Two thirds of that for unbalanced shapes.
    Measured in Nov 2026 on Cascade Lake, Comet Lake, Ice Lake, Meteor Lake,
    Zen 4 and a Sapphire Rapids VM.
*/
static inline int
_nmod_poly_want_karatsuba_one_word(slong len1, slong len2, slong m,
                                   int kind, slong bits, ulong n)
{
    int mid = (kind == _NMOD_POLY_KIND_MID), unknown = 0;
    slong s = mid ? FLINT_MIN(len2, m) : len2;
    slong l = mid ? FLINT_MAX(len2, m)
                  : (kind == _NMOD_POLY_KIND_LOW) ? m : len1;
    int unbalanced = (l >= 2 * s);
    int sizes = (2 * bits + FLINT_BIT_COUNT(len2) <= 50), direct = 0;
    slong lo, hi;

    if (kind == _NMOD_POLY_KIND_SQR)
        lo = (bits < 4) ? WORD_MAX : (bits < 12) ? 13 : 12;
    else
        lo = (bits < 4 && (kind == _NMOD_POLY_KIND_LOW
                           || (kind == _NMOD_POLY_KIND_MUL && unbalanced)))
            ? WORD_MAX : 12;

    if (s < lo)
        return 0;

    if (bits < 14 && !_nmod_poly_karatsuba_small_moduli_ok(s, kind, bits,
                                                           unbalanced))
        return 0;

    if ((kind == _NMOD_POLY_KIND_SQR || kind == _NMOD_POLY_KIND_LOW)
            && s > 16 && s < 20)
        return 0;

#if NMOD_POLY_HAVE_IFMA
    if (len2 <= FLINT_MAX(FLINT_NMOD_POLY_IFMA_MAX_LEN(bits),
                          FLINT_NMOD_POLY_IFMA_BALANCED_MAX_LEN(bits))
            && !_nmod_poly_karatsuba_blocks_ok(s,
                    (kind == _NMOD_POLY_KIND_SQR)
                    ? _NMOD_POLY_KIND_MUL : kind,
                    _NMOD_POLY_KARATSUBA_PAD_STRICT))
        return 0;
#endif

    if (!sizes && bits >= 20 && len1 >= 100)
    {
        /* with n unknown, the smaller of the two bounds */
        unknown = (n == 0);
        direct = unknown || _nmod_poly_fft_small_direct(len1, len2, n);
    }

    if (sizes)
        hi = _nmod_poly_karatsuba_one_prime_max(kind, bits, s, sizes);
    else if (direct)
        hi = 512;
    else
    {
#if NMOD_POLY_HAVE_IFMA
        hi = (kind == _NMOD_POLY_KIND_MUL) ? 2048 : 1024;
#else
        hi = FLINT_NMOD_POLY_KARA_TWO_PRIMES_MAX_LEN;
        if (kind == _NMOD_POLY_KIND_SQR)
            hi /= 2;
#endif
    }

    if (unbalanced)
        hi = 2 * hi / 3;

    return (s <= hi) ? 1 : unknown ? -1 : 0;
}

/*
    The dispatchers (_nmod_poly_mul, _nmod_poly_mullow, _nmod_poly_mulmid)
    work in three steps. Each step calls the next one, an out-of-line function
    as its last action (so that the compiler can make the call a jump, a tail
    call: no new stack frame, and the arguments stay where they are).
    Short products thus pay for one test:
    1. operands shorter than _NMOD_POLY_DISPATCH_MIN_LEN(bits) (and few
       outputs) go to the other algorithms (the functions _other), and so
       do all products in a build without kernels (_NMOD_POLY_HAVE_KERNELS);
    2. _nmod_poly_dispatch decides from the lengths and the bit size of n:
       the other algorithms, the IFMA kernel, the Karatsuba, or "n decides"
       (the functions _kernels);
    3. then _nmod_poly_dispatch_n decides from n: the parity above 52 bits,
       the primes of fft_small, whether one word holds, and the rule of the
       doubles, which changes at n = 1515531528 (the functions _kernels_n).
    The tests check that steps 2 and 3 agree where step 2 decides. In step
    3 the IFMA kernel is called through _nmod_poly_mulmid_ifma_ptr, with mod
    by pointer: with mod by value, gcc 14 and 15 store n and ninv separately
    before jumping to the kernel, which reloads them with one 16-byte load,
    a store-forwarding stall of 5-15 ns. With plain calls instead of these
    jumps, products of 4 to 9 coefficients were 3-10% slower (gcc 13 to 15,
    Nov 2026).
*/
#define _NMOD_POLY_DISPATCH_MIN_LEN(bits) \
    (NMOD_POLY_HAVE_IFMA ? _NMOD_POLY_IFMA_MIN_LEN(bits) : 12)

/* without the IFMA kernel and the Karatsuba (no AVX2), every product goes
   to the other algorithms at once */
#define _NMOD_POLY_HAVE_KERNELS \
    (NMOD_POLY_HAVE_IFMA || NMOD_POLY_HAVE_KARATSUBA)

#define _NMOD_POLY_DISPATCH_OTHER 0
#define _NMOD_POLY_DISPATCH_IFMA 1
#define _NMOD_POLY_DISPATCH_KARATSUBA 2
#define _NMOD_POLY_DISPATCH_BY_N 3

/* the one-word Karatsuba from the bit size of n: 1, 0, or -1 if n decides */
static inline int
_nmod_poly_dispatch_one_word(slong len1, slong len2, slong m, int kind,
                             slong bits)
{
    int mid = (kind == _NMOD_POLY_KIND_MID);

    if (!_nmod_poly_karatsuba_in_range(len2, m, kind, bits))
        return 0;

    if (!_nmod_poly_karatsuba_one_word_bits(bits, mid ? FLINT_MIN(len2, m)
                                                      : len2, len2, mid))
        return -1;

    return _nmod_poly_want_karatsuba_one_word(len1, len2, m, kind, bits, 0);
}

static inline int
_nmod_poly_dispatch(slong len1, slong len2, slong m, int kind, slong bits)
{
    int r;

#if NMOD_POLY_HAVE_IFMA
    int low = (kind == _NMOD_POLY_KIND_LOW);
    int mid = (kind == _NMOD_POLY_KIND_MID);

    if (_nmod_poly_ifma_in_range(len2, m, bits, low))
    {
        if (bits > 52)
            return _nmod_poly_ifma_above_52_in_range(len2, m, bits, low, mid)
                ? _NMOD_POLY_DISPATCH_BY_N : _NMOD_POLY_DISPATCH_OTHER;

        if (bits == 50 && len2 > (low
                ? FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME
                : FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME))
            return _NMOD_POLY_DISPATCH_BY_N;

        /* the one-word Karatsuba first: one block at once, longer blocks
           if its rule takes them */
        if (_nmod_poly_karatsuba_one_block_bits(len1, len2, m, kind, bits))
            return _NMOD_POLY_DISPATCH_KARATSUBA;

        if (_nmod_poly_karatsuba_before_ifma_bits(len1, len2, m, kind, bits))
        {
            r = _nmod_poly_dispatch_one_word(len1, len2, m, kind, bits);

            if (r != 0)
                return (r > 0) ? _NMOD_POLY_DISPATCH_KARATSUBA
                               : _NMOD_POLY_DISPATCH_BY_N;
        }

        /* beyond FLINT_NMOD_POLY_IFMA_MAX_LEN (balanced shapes), n decides
           (_nmod_poly_want_ifma) */
        return (!low && m >= len2 && len2 > FLINT_NMOD_POLY_IFMA_MAX_LEN(bits))
            ? _NMOD_POLY_DISPATCH_BY_N : _NMOD_POLY_DISPATCH_IFMA;
    }
#endif

    if (bits > 30)
        return (NMOD_POLY_HAVE_KARATSUBA_FP50
                && _nmod_poly_karatsuba_in_range(len2, m, kind, bits))
            ? _NMOD_POLY_DISPATCH_BY_N : _NMOD_POLY_DISPATCH_OTHER;

    r = _nmod_poly_dispatch_one_word(len1, len2, m, kind, bits);

    return (r > 0) ? _NMOD_POLY_DISPATCH_KARATSUBA
         : (r == 0) ? _NMOD_POLY_DISPATCH_OTHER : _NMOD_POLY_DISPATCH_BY_N;
}

/* where n decides: the IFMA kernel, with the one-word Karatsuba before it
   as in _nmod_poly_dispatch, else the Karatsuba if its rule takes the
   product, else the other algorithms */
static inline int
_nmod_poly_dispatch_n(slong len1, slong len2, slong m, int kind, nmod_t mod)
{
    int ifma = _nmod_poly_want_ifma(len1, len2, m, mod,
                                    kind == _NMOD_POLY_KIND_LOW);

    if (ifma && _nmod_poly_karatsuba_one_block_bits(len1, len2, m, kind,
                                                    NMOD_BITS(mod)))
        return _NMOD_POLY_DISPATCH_KARATSUBA;

    if ((!ifma || _nmod_poly_karatsuba_before_ifma(len1, len2, m, kind, mod))
            && _nmod_poly_karatsuba_in_range(len2, m, kind, NMOD_BITS(mod))
            && _nmod_poly_want_karatsuba(len1, len2, m, kind, mod.n))
        return _NMOD_POLY_DISPATCH_KARATSUBA;

    return ifma ? _NMOD_POLY_DISPATCH_IFMA : _NMOD_POLY_DISPATCH_OTHER;
}

void _nmod_poly_inv_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
void _nmod_poly_div_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr P, slong Plen, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
int nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n);
void _nmod_poly_divrem_q0_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, nn_srcptr B, slong lenA, ulong invL, nmod_t mod);

/* product with a precomputed transform, see geometric_progression.c */
struct mul_precomp_struct;
int _nmod_geometric_progression_mulmid_fft(nn_ptr z, slong zl, slong zh,
        nn_srcptr a, slong an, const struct mul_precomp_struct * F,
        slong s, slong bn, nmod_t mod);

#endif
