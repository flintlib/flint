/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/*
    Transform-space versions of the matrix operations of the half-gcd.

    The 2 x 2 products have four outputs, each a sum of two products, so
    with fft_small every operand is transformed once and every output is
    transformed back once:

      matrix_mul (M <- M M1)        8 forward, 4 inverse transforms
                                    (Strassen-like: 7 products, 21)
      matrix_adjust (low parts)     6 forward, 2 inverse (4 products: 12)
      matrix_adjust_mul (both of    10 forward, 6 inverse, the entries of
        the above, sharing M1)      M1 transformed once (separately: 20)
      mul_matrix_vector (gcdext)    6 forward, 2 inverse (4 products: 12)
      matrix_apply (mod B^N - 1)    7 forward, 2 inverse cyclic transforms
                                    of half length (4 cyclic products: 12)

    The entries are nonnegative. Outputs that are sums of two products are
    exported unsigned; the differences of matrix_adjust are exported
    signed; the differences of matrix_apply are made nonnegative modulo
    B^N - 1 by replacing -x with (B^N - 1) - x, whose transform is the
    transform of B^N - 1 (all ones) minus that of x.
*/

/* Whether products of lengths at least an >= bn should go through
   fft_small when the transforms are shared as above: the rule of
   flint_mpn_mul with the threshold T replaced by the given cutoff. */
static int
use_fft(mp_size_t an, mp_size_t bn, mp_size_t cutoff)
{
    if (an < bn)
        FLINT_SWAP(mp_size_t, an, bn);
    return (bn >= cutoff) || (an + bn > 2 * cutoff && 3 * bn > cutoff);
}

#ifndef FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF
# define FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF (3 * FLINT_FFT_SMALL_MUL_THRESHOLD / 10)
#endif
#ifndef FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF
# define FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF (3 * FLINT_FFT_SMALL_MUL_THRESHOLD / 10)
#endif
#ifndef FLINT_MPN_HGCD_APPLY_FFT_CUTOFF
# define FLINT_MPN_HGCD_APPLY_FFT_CUTOFF (3 * FLINT_FFT_SMALL_MUL_THRESHOLD / 10)
#endif

/* Whether count operands of the plan fit within the global bound on
   transform storage; if not, the functions below decline, and their
   callers fall back to separate products, which need less memory. */
static int
pool_fits(const fft_small_plan_struct * P, ulong count)
{
    return count * fft_small_op_sizeof_data(P) <= flint_fft_small_max_transformed_ring_size;
}

/* z (zn limbs) = the exact nonnegative value of Z, a sum of at most terms
   products of operands whose products have at most nslots slots; the
   value must fit in zn limbs. */
static void
export_unsigned(mp_ptr z, mp_size_t zn, fft_small_op_struct * Z,
    ulong terms, ulong nslots, const fft_small_plan_struct * P)
{
    ulong need = fft_small_export_mpn_limbs(P, terms, nslots);
    ulong en = FLINT_MIN((ulong) zn, need);

    fft_small_ifft_trunc(Z, nslots, n_cdiv(en * FLINT_BITS, P->bits), NULL, P);
    fft_small_export_mpn(z, en, Z, P);
    if ((ulong) zn > en)
        flint_mpn_zero(z + en, zn - en);
}

/* (buf, returned length) = |Z|, sign in *sgn; buf needs
   fft_small_export_mpn_signed_limbs(P, nslots) limbs. */
static mp_size_t
export_signed(mp_ptr buf, int * sgn, fft_small_op_struct * Z,
    ulong nslots, const fft_small_plan_struct * P)
{
    ulong need = fft_small_export_mpn_signed_limbs(P, nslots);
    mp_size_t len;

    fft_small_ifft_trunc(Z, nslots, 0, NULL, P);
    fft_small_export_mpn_signed(buf, need, sgn, Z, nslots, P);

    len = need;
    MPN_NORM(buf, len);
    return len;
}

/* The final stage of M <- M M1 with the entries of M of rn limbs and
   those of M1 of mn limbs: the new entries (r00, r01, r10, r11) are
   exported from Z00, Z01, Z10, Z11 and M->n is set. */
static void
matmul_finish(flint_mpn_hgcd_matrix_struct * M,
    fft_small_op_struct * Z00, fft_small_op_struct * Z01,
    fft_small_op_struct * Z10, fft_small_op_struct * Z11,
    ulong nslots, mp_size_t rn, mp_size_t mn, const fft_small_plan_struct * P)
{
    /* The entries of the product have at most rn + mn + 1 limbs; the
       allocation of M allows this. */
    mp_size_t zn = rn + mn + 1, n;
    FLINT_ASSERT(zn <= M->alloc);

    export_unsigned(M->p[0][0], zn, Z00, 2, nslots, P);
    export_unsigned(M->p[0][1], zn, Z01, 2, nslots, P);
    export_unsigned(M->p[1][0], zn, Z10, 2, nslots, P);
    export_unsigned(M->p[1][1], zn, Z11, 2, nslots, P);

    /* Index of last potentially non-zero limb, size is one greater. For
       consecutive half-gcd matrices this loop runs at most three times
       (see _flint_mpn_hgcd_matrix_mul), but other products may shrink
       more. */
    n = rn + mn;
    while (n > 0 && (M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) == 0)
        n--;

    FLINT_ASSERT((M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) > 0);

    M->n = n + 1;
}

/* (a; b) <- a_hi B^p + S_a, b_hi B^p + S_b, the final stage of the matrix
   adjustment: S_a, S_b (signed, of at most p + mn limbs) are exported
   from Za and Zb, whose products have at most nslots slots. a and b have
   n limbs (and room for n + 1). Returns the new size. buf needs
   fft_small_export_mpn_signed_limbs(P, nslots) limbs. */
static mp_size_t
adjust_finish(fft_small_op_struct * Za, fft_small_op_struct * Zb,
    ulong nslots, mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p,
    mp_size_t mn, mp_ptr buf, const fft_small_plan_struct * P)
{
    mp_limb_t ah, bh;
    mp_size_t sn;
    int sgn;

    /* a = a_hi B^p + S */
    sn = export_signed(buf, &sgn, Za, nslots, P);
    FLINT_ASSERT(sn <= p + mn);
    flint_mpn_zero(ap, p);
    if (sn == 0)
        ah = 0;
    else if (!sgn)
        ah = mpn_add(ap, ap, n, buf, sn);
    else
    {
        mp_limb_t cy = mpn_sub(ap, ap, n, buf, sn);
        FLINT_ASSERT(cy == 0);
        (void) cy;
        ah = 0;
    }

    sn = export_signed(buf, &sgn, Zb, nslots, P);
    FLINT_ASSERT(sn <= p + mn);
    flint_mpn_zero(bp, p);
    if (sn == 0)
        bh = 0;
    else if (!sgn)
        bh = mpn_add(bp, bp, n, buf, sn);
    else
    {
        mp_limb_t cy = mpn_sub(bp, bp, n, buf, sn);
        FLINT_ASSERT(cy == 0);
        (void) cy;
        bh = 0;
    }

    (void) mn;

    if (ah > 0 || bh > 0)
    {
        ap[n] = ah;
        bp[n] = bh;
        n++;
    }
    else
    {
        /* The subtraction can reduce the size by at most one limb. */
        if (ap[n - 1] == 0 && bp[n - 1] == 0)
            n--;
    }
    FLINT_ASSERT(ap[n - 1] > 0 || bp[n - 1] > 0);
    return n;
}

/******************************************************************************

    M <- M * M1

******************************************************************************/

int
_flint_mpn_hgcd_matrix_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], Y[2][2], T;
    mp_size_t rn = M->n, mn = M1->n;
    ulong nslots;
    int i, j;

    if (!use_fft(rn, mn, FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF))
        return 0;

    R = get_default_mpn_ctx();

    /* Each output is a sum of two products of nonnegative entries. */
    if (!fft_small_plan_init_mpn(P, R, rn, mn, 2, 0))
        return 0;

    if (!pool_fits(P, 9))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, rn, mn);

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 9, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_op_init_borrowed(Y[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_fft_mpn(X[i] + j, M->p[i][j], rn, P);
            fft_small_fft_mpn(Y[i] + j, M1->p[i][j], mn, P);
        }

    fft_small_op_init_borrowed(&T, P, fft_small_op_pool_take(pool));

    /*
        r00' = r00 m00 + r01 m10  -> T
        r01' = r00 m01 + r01 m11  -> X00 (r00, r01 dead)
        r10' = r10 m00 + r11 m10  -> X01
        r11' = r10 m01 + r11 m11  -> X10 (r10, r11 dead)
    */
    fft_small_op_mul(&T, X[0] + 0, Y[0] + 0, P);
    fft_small_op_addmul(&T, X[0] + 1, Y[1] + 0, P);
    fft_small_op_mul(X[0] + 0, X[0] + 0, Y[0] + 1, P);
    fft_small_op_addmul(X[0] + 0, X[0] + 1, Y[1] + 1, P);
    fft_small_op_mul(X[0] + 1, X[1] + 0, Y[0] + 0, P);
    fft_small_op_addmul(X[0] + 1, X[1] + 1, Y[1] + 0, P);
    fft_small_op_mul(X[1] + 0, X[1] + 0, Y[0] + 1, P);
    fft_small_op_addmul(X[1] + 0, X[1] + 1, Y[1] + 1, P);

    matmul_finish(M, &T, X[0] + 0, X[0] + 1, X[1] + 0, nslots, rn, mn, P);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_pool_give(pool, X[i][j].data);
            fft_small_op_pool_give(pool, Y[i][j].data);
        }
    fft_small_op_pool_give(pool, T.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);
    return 1;
}

/******************************************************************************

    (a; b) <- M^{-1} (a; b) on the low p limbs

******************************************************************************/

/*
    Sets a = a_hi B^p + (m11 a_lo - m01 b_lo), b = b_hi B^p + (m00 b_lo -
    m10 a_lo) where a_lo, b_lo are the low p limbs; see
    _flint_mpn_hgcd_matrix_adjust. Returns the new size, or 0 if fft_small
    is not used.
*/
mp_size_t
_flint_mpn_hgcd_matrix_adjust_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B;
    mp_size_t mn = M->n, an, bn;
    mp_ptr buf;
    ulong nslots, need;
    int i, j;

    if (!use_fft(p, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF))
        return 0;

    an = bn = p;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    /* Each output is a difference of two products. */
    if (!fft_small_plan_init_mpn(P, R, p, mn, 2, 1))
        return 0;

    if (!pool_fits(P, 6))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, p, mn);
    need = fft_small_export_mpn_signed_limbs(P, nslots);

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 6, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_fft_mpn(X[i] + j, M->p[i][j], mn, P);
        }

    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_fft_mpn(&A, ap, an, P);
    fft_small_fft_mpn(&B, bp, bn, P);

    /* a' = m11 a - m01 b -> X11, b' = m00 b - m10 a -> X00 */
    fft_small_op_mul(X[1] + 1, X[1] + 1, &A, P);
    fft_small_op_submul(X[1] + 1, X[0] + 1, &B, P);
    fft_small_op_mul(X[0] + 0, X[0] + 0, &B, P);
    fft_small_op_submul(X[0] + 0, X[1] + 0, &A, P);

    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_give(pool, X[0][1].data);
    fft_small_op_pool_give(pool, X[1][0].data);

    buf = flint_malloc(need * sizeof(mp_limb_t));
    n = adjust_finish(X[1] + 1, X[0] + 0, nslots, n, ap, bp, p, mn, buf, P);
    flint_free(buf);

    fft_small_op_pool_give(pool, X[0][0].data);
    fft_small_op_pool_give(pool, X[1][1].data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    return n;
}

/******************************************************************************

    (a; b) <- M1^{-1} (a; b) on the low p limbs, then M <- M * M1

******************************************************************************/

/*
    The two operations following the second recursive call of the half-gcd,
    which share the matrix M1: with entries of M of about the size of a_lo
    and b_lo, one plan serves both, and the entries of M1 are transformed
    once (10 forward transforms instead of 14, out of 20 in total).
    Returns the new size of a and b, or 0 (doing nothing) if fft_small is
    not used, in which case the caller performs the two operations
    separately.
*/
mp_size_t
_flint_mpn_hgcd_matrix_adjust_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], Y[2][2], A, B, Za, Zb;
    mp_size_t rn = M->n, mn = M1->n, an, bn, lmax, lmin;
    mp_ptr buf;
    ulong nslots_adj, nslots_mul, need;
    int i, j;

    if (!use_fft(p, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF) ||
        !use_fft(rn, mn, FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF))
        return 0;

    /* one transform length for both products: only when it does not
       penalize the shorter one much */
    lmax = FLINT_MAX(p, rn) + mn;
    lmin = FLINT_MIN(p, rn) + mn;
    if (4 * lmax > 5 * lmin)
        return 0;

    an = bn = p;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    /* differences of two products for the adjustment, sums of two products
       for the matrix product */
    if (!fft_small_plan_init_mpn(P, R, FLINT_MAX(p, rn), mn, 2, 1))
        return 0;

    if (!pool_fits(P, 12))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots_adj = fft_small_plan_mpn_nslots(P, p, mn);
    nslots_mul = fft_small_plan_mpn_nslots(P, rn, mn);
    need = fft_small_export_mpn_signed_limbs(P, nslots_adj);

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 12, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_op_init_borrowed(Y[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_fft_mpn(Y[i] + j, M1->p[i][j], mn, P);
            fft_small_fft_mpn(X[i] + j, M->p[i][j], rn, P);
        }

    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&Za, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&Zb, P, fft_small_op_pool_take(pool));
    fft_small_fft_mpn(&A, ap, an, P);
    fft_small_fft_mpn(&B, bp, bn, P);

    /* a' = m11 a - m01 b -> Za, b' = m00 b - m10 a -> Zb */
    fft_small_op_mul(&Za, Y[1] + 1, &A, P);
    fft_small_op_submul(&Za, Y[0] + 1, &B, P);
    fft_small_op_mul(&Zb, Y[0] + 0, &B, P);
    fft_small_op_submul(&Zb, Y[1] + 0, &A, P);

    /*
        r00' = r00 m00 + r01 m10  -> A (dead)
        r01' = r00 m01 + r01 m11  -> X00 (r00, r01 dead)
        r10' = r10 m00 + r11 m10  -> X01
        r11' = r10 m01 + r11 m11  -> X10 (r10, r11 dead)
    */
    fft_small_op_mul(&A, X[0] + 0, Y[0] + 0, P);
    fft_small_op_addmul(&A, X[0] + 1, Y[1] + 0, P);
    fft_small_op_mul(X[0] + 0, X[0] + 0, Y[0] + 1, P);
    fft_small_op_addmul(X[0] + 0, X[0] + 1, Y[1] + 1, P);
    fft_small_op_mul(X[0] + 1, X[1] + 0, Y[0] + 0, P);
    fft_small_op_addmul(X[0] + 1, X[1] + 1, Y[1] + 0, P);
    fft_small_op_mul(X[1] + 0, X[1] + 0, Y[0] + 1, P);
    fft_small_op_addmul(X[1] + 0, X[1] + 1, Y[1] + 1, P);

    buf = flint_malloc(need * sizeof(mp_limb_t));
    n = adjust_finish(&Za, &Zb, nslots_adj, n, ap, bp, p, mn, buf, P);
    flint_free(buf);

    matmul_finish(M, &A, X[0] + 0, X[0] + 1, X[1] + 0, nslots_mul, rn, mn, P);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_pool_give(pool, X[i][j].data);
            fft_small_op_pool_give(pool, Y[i][j].data);
        }
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_give(pool, Za.data);
    fft_small_op_pool_give(pool, Zb.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    return n;
}

/******************************************************************************

    (r; b) <- (a; b) M   (cofactor update in gcdext)

******************************************************************************/

/*
    Sets (r, b) = (m00 a + m10 b, m01 a + m11 b) where a, b have n limbs
    (not necessarily normalized); r and b need room for n + M->n + 1 limbs.
    Returns the new size (normalized), or 0 if fft_small is not used.
*/
mp_size_t
_flint_mpn_hgcd_mul_matrix_vector_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B;
    mp_size_t mn = M->n, an, bn, zn;
    ulong nslots;
    int i, j;

    if (!use_fft(n, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF))
        return 0;

    an = bn = n;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    if (!fft_small_plan_init_mpn(P, R, n, mn, 2, 0))
        return 0;

    if (!pool_fits(P, 6))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, n, mn);

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 6, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_fft_mpn(X[i] + j, M->p[i][j], mn, P);
        }

    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_fft_mpn(&A, ap, an, P);
    fft_small_fft_mpn(&B, bp, bn, P);

    /* r = m00 a + m10 b -> X00, b = m01 a + m11 b -> X01 */
    fft_small_op_mul(X[0] + 0, X[0] + 0, &A, P);
    fft_small_op_addmul(X[0] + 0, X[1] + 0, &B, P);
    fft_small_op_mul(X[0] + 1, X[0] + 1, &A, P);
    fft_small_op_addmul(X[0] + 1, X[1] + 1, &B, P);

    /* The products have at most n + mn + 1 limbs, but like the classical
       code we may only write limb n + mn when it is nonzero (the caller's
       allocation is sized for the normalized result), so the exports go
       through a temporary. */
    zn = n + mn + 1;
    {
        mp_ptr t = flint_malloc(2 * zn * sizeof(mp_limb_t));

        export_unsigned(t, zn, X[0] + 0, 2, nslots, P);
        export_unsigned(t + zn, zn, X[0] + 1, 2, nslots, P);

        flint_mpn_copyi(rp, t, zn - 1);
        flint_mpn_copyi(bp, t + zn, zn - 1);
        if ((t[zn - 1] | t[2 * zn - 1]) != 0)
        {
            rp[zn - 1] = t[zn - 1];
            bp[zn - 1] = t[2 * zn - 1];
        }
        else
            zn--;

        flint_free(t);
    }

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_pool_give(pool, X[i][j].data);
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    while ((rp[zn - 1] | bp[zn - 1]) == 0)
        zn--;

    return zn;
}

/******************************************************************************

    (a; b) <- M^{-1} (a; b), result known to fit in nn limbs

******************************************************************************/

/* r = v mod (B^N - 1) in [0, B^N - 1), for v of vn limbs; r may equal v */
static void
fold_bnm1(mp_ptr r, mp_srcptr v, mp_size_t vn, mp_size_t N)
{
    mp_limb_t cy;
    mp_size_t m, i;

    if (vn <= N)
    {
        if (r != v)
            flint_mpn_copyi(r, v, vn);
        flint_mpn_zero(r + vn, N - vn);
    }
    else
    {
        if (r != v)
            flint_mpn_copyi(r, v, N);
        v += N;
        vn -= N;

        while (vn > 0)
        {
            m = FLINT_MIN(vn, N);
            cy = mpn_add(r, r, N, v, m);
            while (cy)
                cy = mpn_add_1(r, r, N, cy);   /* B^N = 1 */
            v += m;
            vn -= m;
        }
    }

    for (i = 0; i < N; i++)
        if (r[i] != UWORD_MAX)
            return;
    flint_mpn_zero(r, N);                      /* B^N - 1 = 0 */
}

int
_flint_mpn_hgcd_matrix_apply_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t mn[2][2], mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t nn)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B, T;
    ulong N, zn;
    mp_ptr t, z;
    int i, j;

    if (!use_fft(nn, M->n, FLINT_MPN_HGCD_APPLY_FFT_CUTOFF))
        return 0;

    R = get_default_mpn_ctx();

    /* The results are nonnegative and less than B^nn, so they are
       determined by their residues mod B^N - 1 for N > nn. Each output
       is a sum of two cyclic products. */
    N = nn + 1;
    if (!fft_small_plan_init_mpn_cyclic(P, R, &N, 2))
        return 0;

    if (!pool_fits(P, 7))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    FLINT_ASSERT(N > (ulong) nn);

    /* the matrix entries must fit (they do in hgcd_reduce), as must a and
       b after one folding */
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            FLINT_ASSERT(mn[i][j] > 0);
    if ((ulong) FLINT_MAX(FLINT_MAX(mn[0][0], mn[0][1]),
                          FLINT_MAX(mn[1][0], mn[1][1])) > N)
    {
        fft_small_plan_clear(P);
        return 0;
    }

    zn = N + (2 * P->bits + P->depth + FLINT_BIT_COUNT(2)) / FLINT_BITS + 1;

    t = flint_malloc((N + zn) * sizeof(mp_limb_t));
    z = t + N;

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 7, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_fft_mpn(X[i] + j, M->p[i][j], mn[i][j], P);
        }

    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&T, P, fft_small_op_pool_take(pool));

    /* a, b mod B^N - 1 */
    fold_bnm1(t, ap, n, N);
    fft_small_fft_mpn(&A, t, N, P);
    fold_bnm1(t, bp, n, N);
    fft_small_fft_mpn(&B, t, N, P);

    /* transform of B^N - 1 */
    for (i = 0; i < (slong) N; i++)
        t[i] = UWORD_MAX;
    fft_small_fft_mpn(&T, t, N, P);

    /*
        a' = m11 a + m01 ((B^N - 1) - b) -> X11
        b' = m00 b + m10 ((B^N - 1) - a) -> X00
    */
    fft_small_op_mul(X[1] + 1, X[1] + 1, &A, P);
    fft_small_op_mul(X[0] + 0, X[0] + 0, &B, P);
    fft_small_op_sub(&A, &T, &A, P);     /* (B^N - 1) - a */
    fft_small_op_sub(&B, &T, &B, P);     /* (B^N - 1) - b */
    fft_small_op_addmul(X[1] + 1, X[0] + 1, &B, P);
    fft_small_op_addmul(X[0] + 0, X[1] + 0, &A, P);

    fft_small_ifft(X[1] + 1, P);
    fft_small_export_mpn(z, zn, X[1] + 1, P);
    fold_bnm1(t, z, zn, N);
    FLINT_ASSERT(flint_mpn_zero_p(t + nn, N - nn));
    flint_mpn_copyi(ap, t, nn);

    fft_small_ifft(X[0] + 0, P);
    fft_small_export_mpn(z, zn, X[0] + 0, P);
    fold_bnm1(t, z, zn, N);
    FLINT_ASSERT(flint_mpn_zero_p(t + nn, N - nn));
    flint_mpn_copyi(bp, t, nn);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_pool_give(pool, X[i][j].data);
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_give(pool, T.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);
    flint_free(t);

    return 1;
}

#endif
