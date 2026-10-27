/*
    From GMP's mpn/generic/hgcd.c, hgcd_appr.c, hgcd_reduce.c, hgcd_step.c:
    Copyright 2003-2005, 2008, 2011, 2012 Free Software Foundation, Inc.

    From GMP's mpn/generic/hgcd_matrix.c, matrix22_mul.c:
    Copyright 2003-2005, 2008, 2012 Free Software Foundation, Inc.

    From GMP's mpn/generic/gcd_subdiv_step.c:
    Copyright 2003-2005, 2008, 2010, 2011 Free Software Foundation, Inc.

    Contributed to the GNU project by Niels Möller and Marco Bodrato.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

/*
    Half-gcd (HGCD) following N. Möller, On Schönhage's algorithm and
    subquadratic integer gcd computation, Math. Comp. 77 (2008), 589-607,
    ported from GMP. The algorithms and the storage analysis are GMP's;
    FLINT's multiplication and division are used throughout, and the
    matrix operations dispatch to fft_small transform-space versions
    (hgcd_fft.c) for large operands.
*/

#define MUL(rp, ap, an, bp, bn) \
    do { \
        if ((an) >= (bn)) \
            flint_mpn_mul(rp, ap, an, bp, bn); \
        else \
            flint_mpn_mul(rp, bp, bn, ap, an); \
    } while (0)

/* Extract a full limb from (xh, xl) after a left shift by count bits;
   valid for 0 <= count < FLINT_BITS. */
#define EXTRACT_LIMB(count, xh, xl) \
    (((xh) << (count)) | (((xl) >> 1) >> ((FLINT_BITS - 1 - (count)))))

/******************************************************************************

    Matrix operations

******************************************************************************/

/* For input of size n, matrix elements are of size at most ceil(n/2) - 1,
   but we need two limbs extra. */
void
_flint_mpn_hgcd_matrix_init(flint_mpn_hgcd_matrix_struct * M, mp_size_t n, mp_ptr p)
{
    mp_size_t s = (n + 1) / 2 + 1;
    M->alloc = s;
    M->n = 1;
    flint_mpn_zero(p, 4 * s);
    M->p[0][0] = p;
    M->p[0][1] = p + s;
    M->p[1][0] = p + 2 * s;
    M->p[1][1] = p + 3 * s;

    M->p[0][0][0] = M->p[1][1][0] = 1;
}

/* Update column col, adding in Q * column (1 - col). Temporary storage:
   qn + n <= M->alloc, where n is the size of the largest element in column
   1 - col. */
void
_flint_mpn_hgcd_matrix_update_q(flint_mpn_hgcd_matrix_struct * M,
    mp_srcptr qp, mp_size_t qn, unsigned int col, mp_ptr tp)
{
    FLINT_ASSERT(col < 2);

    if (qn == 1)
    {
        mp_limb_t q = qp[0];
        mp_limb_t c0, c1;

        c0 = mpn_addmul_1(M->p[0][col], M->p[0][1 - col], M->n, q);
        c1 = mpn_addmul_1(M->p[1][col], M->p[1][1 - col], M->n, q);

        M->p[0][col][M->n] = c0;
        M->p[1][col][M->n] = c1;

        M->n += (c0 | c1) != 0;
    }
    else
    {
        unsigned int row;

        /* Carries for the unlikely case that we get both high words from
           the multiplication and carries from the addition. */
        mp_limb_t c[2];
        mp_size_t n;

        /* The matrix will not necessarily grow in size by qn, so we need
           normalization in order not to overflow M. */
        for (n = M->n; n + qn > M->n; n--)
        {
            FLINT_ASSERT(n > 0);
            if (M->p[0][1 - col][n - 1] > 0 || M->p[1][1 - col][n - 1] > 0)
                break;
        }

        FLINT_ASSERT(qn + n <= M->alloc);

        for (row = 0; row < 2; row++)
        {
            MUL(tp, M->p[row][1 - col], n, qp, qn);

            FLINT_ASSERT(n + qn >= M->n);
            c[row] = mpn_add(M->p[row][col], tp, n + qn, M->p[row][col], M->n);
        }

        n += qn;

        if (c[0] | c[1])
        {
            M->p[0][col][n] = c[0];
            M->p[1][col][n] = c[1];
            n++;
        }
        else
        {
            n -= (M->p[0][col][n - 1] | M->p[1][col][n - 1]) == 0;
            FLINT_ASSERT(n >= M->n);
        }
        M->n = n;
    }

    FLINT_ASSERT(M->n < M->alloc);
}

/* Multiply M by M1 from the right. Since the M1 elements fit in
   FLINT_BITS - 1 bits, M grows by at most one limb. The temporary space
   tp (M->n limbs in GMP) is not used. */
void
_flint_mpn_hgcd_matrix_mul_1(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix1_struct * M1, mp_ptr tp)
{
    mp_size_t n0, n1;

    (void) tp;

    /* in place: the first output may overwrite the first input */
    n0 = _flint_mpn_hgcd_mul_matrix1_vector(M1, M->p[0][0], M->p[0][0], M->p[0][1], M->n);
    n1 = _flint_mpn_hgcd_mul_matrix1_vector(M1, M->p[1][0], M->p[1][0], M->p[1][1], M->n);

    /* Depends on zero initialization */
    M->n = FLINT_MAX(n0, n1);
    FLINT_ASSERT(M->n < M->alloc);
}

/* 2 x 2 matrix products (port of GMP's mpn_matrix22_mul) */

/* Inputs are unsigned. */
static int
abs_sub_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
    int c = mpn_cmp(ap, bp, n);
    if (c >= 0)
    {
        mpn_sub_n(rp, ap, bp, n);
        return 0;
    }
    else
    {
        mpn_sub_n(rp, bp, ap, n);
        return 1;
    }
}

static int
add_signed_n(mp_ptr rp, mp_srcptr ap, int as, mp_srcptr bp, int bs, mp_size_t n)
{
    if (as != bs)
        return as ^ abs_sub_n(rp, ap, bp, n);
    else
    {
        mp_limb_t cy = mpn_add_n(rp, ap, bp, n);
        FLINT_ASSERT(cy == 0);
        (void) cy;
        return as;
    }
}

mp_size_t
_flint_mpn_hgcd_matrix22_mul_strassen_itch(mp_size_t rn, mp_size_t mn)
{
    if (rn < FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF
        || mn < FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF)
        return 3 * rn + 2 * mn;
    else
        return 3 * (rn + mn) + 5;
}

/*
    Algorithm (Bodrato, "A Strassen-like Matrix Multiplication suited for
    squaring and higher power computation", ISSAC 2010):

    / s0 \   /  1  0  0  0 \ / r0 \
    | s1 |   |  0  1  0  1 | | r1 |
    | s2 |   |  0  0 -1  1 | | r2 |
    | s3 | = |  0  1 -1  1 | \ r3 /
    | s4 |   | -1  1 -1  1 |
    | s5 |   |  0  1  0  0 |
    \ s6 /   \  0  0  1  0 /

    / t0 \   /  1  0  0  0 \ / m0 \
    | t1 |   |  0  1  0  1 | | m1 |
    | t2 |   |  0  0 -1  1 | | m2 |
    | t3 | = |  0  1 -1  1 | \ m3 /
    | t4 |   | -1  1 -1  1 |
    | t5 |   |  0  1  0  0 |
    \ t6 /   \  0  0  1  0 /

    / r0 \   / 1 0  0  0  0  1  0 \ / s0*t0 \
    | r1 | = | 0 0 -1  1 -1  1  0 | | s1*t1 |
    | r2 |   | 0 1  0 -1  0 -1 -1 | | s2*t2 |
    \ r3 /   \ 0 1  1 -1  0 -1  0 / | s3*t3 |
                                    | s4*t5 |
                                    | s5*t6 |
                                    \ s6*t4 /

    Computes R = R * M. Resulting elements are of size up to rn + mn + 1.
    Temporary storage: 3 rn + 3 mn + 5.
*/
static void
matrix22_mul_strassen(mp_ptr r0, mp_ptr r1, mp_ptr r2, mp_ptr r3, mp_size_t rn,
    mp_srcptr m0, mp_srcptr m1, mp_srcptr m2, mp_srcptr m3, mp_size_t mn,
    mp_ptr tp)
{
    mp_ptr s0, t0, u0, u1;
    int r1s, r3s, s0s, t0s, u1s;
    s0 = tp; tp += rn + 1;
    t0 = tp; tp += mn + 1;
    u0 = tp; tp += rn + mn + 1;
    u1 = tp; /* rn + mn + 2 */

    MUL(u0, r1, rn, m2, mn);            /* u5 = s5 * t6 */
    r3s = abs_sub_n(r3, r3, r2, rn);    /* r3 - r2 */
    if (r3s)
    {
        r1s = abs_sub_n(r1, r1, r3, rn);
        r1[rn] = 0;
    }
    else
    {
        r1[rn] = mpn_add_n(r1, r1, r3, rn);
        r1s = 0;                        /* r1 - r2 + r3  */
    }
    if (r1s)
    {
        s0[rn] = mpn_add_n(s0, r1, r0, rn);
        s0s = 0;
    }
    else if (r1[rn] != 0)
    {
        s0[rn] = r1[rn] - mpn_sub_n(s0, r1, r0, rn);
        s0s = 1;                        /* s4 = -r0 + r1 - r2 + r3 */
                                        /* Reverse sign! */
    }
    else
    {
        s0s = abs_sub_n(s0, r0, r1, rn);
        s0[rn] = 0;
    }
    MUL(u1, r0, rn, m0, mn);            /* u0 = s0 * t0 */
    r0[rn + mn] = mpn_add_n(r0, u0, u1, rn + mn);
    FLINT_ASSERT(r0[rn + mn] < 2);      /* u0 + u5 */

    t0s = abs_sub_n(t0, m3, m2, mn);
    u1s = r3s ^ t0s ^ 1;                /* Reverse sign! */
    MUL(u1, r3, rn, t0, mn);            /* u2 = s2 * t2 */
    u1[rn + mn] = 0;
    if (t0s)
    {
        t0s = abs_sub_n(t0, m1, t0, mn);
        t0[mn] = 0;
    }
    else
    {
        t0[mn] = mpn_add_n(t0, t0, m1, mn);
    }

    if (t0[mn] != 0)
    {
        MUL(r3, r1, rn, t0, mn + 1);    /* u3 = s3 * t3 */
        FLINT_ASSERT(r1[rn] < 2);
        if (r1[rn] != 0)
            mpn_add_n(r3 + rn, r3 + rn, t0, mn + 1);
    }
    else
    {
        MUL(r3, r1, rn + 1, t0, mn);
    }

    FLINT_ASSERT(r3[rn + mn] < 4);

    u0[rn + mn] = 0;
    if (r1s ^ t0s)
    {
        r3s = abs_sub_n(r3, u0, r3, rn + mn + 1);
    }
    else
    {
        mpn_add_n(r3, r3, u0, rn + mn + 1);
        r3s = 0;                        /* u3 + u5 */
    }

    if (t0s)
    {
        t0[mn] = mpn_add_n(t0, t0, m0, mn);
    }
    else if (t0[mn] != 0)
    {
        t0[mn] -= mpn_sub_n(t0, t0, m0, mn);
    }
    else
    {
        t0s = abs_sub_n(t0, t0, m0, mn);
    }
    MUL(u0, r2, rn, t0, mn + 1);        /* u6 = s6 * t4 */
    FLINT_ASSERT(u0[rn + mn] < 2);
    if (r1s)
    {
        mpn_sub_n(r1, r2, r1, rn);
    }
    else
    {
        r1[rn] += mpn_add_n(r1, r1, r2, rn);
    }
    rn++;
    t0s = add_signed_n(r2, r3, r3s, u0, t0s, rn + mn);
                                        /* u3 + u5 + u6 */
    FLINT_ASSERT(r2[rn + mn - 1] < 4);
    r3s = add_signed_n(r3, r3, r3s, u1, u1s, rn + mn);
                                        /* -u2 + u3 + u5  */
    FLINT_ASSERT(r3[rn + mn - 1] < 3);
    MUL(u0, s0, rn, m1, mn);            /* u4 = s4 * t5 */
    FLINT_ASSERT(u0[rn + mn - 1] < 2);
    t0[mn] = mpn_add_n(t0, m3, m1, mn);
    MUL(u1, r1, rn, t0, mn + 1);        /* u1 = s1 * t1 */
    mn += rn;
    FLINT_ASSERT(u1[mn - 1] < 4);
    FLINT_ASSERT(u1[mn] == 0);
    add_signed_n(r1, r3, r3s, u0, s0s, mn);
                                        /* -u2 + u3 - u4 + u5  */
    FLINT_ASSERT(r1[mn - 1] < 2);
    if (r3s)
    {
        mpn_add_n(r3, u1, r3, mn);
    }
    else
    {
        mpn_sub_n(r3, u1, r3, mn);
                                        /* u1 + u2 - u3 - u5  */
    }
    FLINT_ASSERT(r3[mn - 1] < 2);
    if (t0s)
    {
        mpn_add_n(r2, u1, r2, mn);
    }
    else
    {
        mpn_sub_n(r2, u1, r2, mn);
                                        /* u1 - u3 - u5 - u6  */
    }
    FLINT_ASSERT(r2[mn - 1] < 2);
}

void
_flint_mpn_hgcd_matrix22_mul_strassen(mp_ptr r0, mp_ptr r1, mp_ptr r2, mp_ptr r3, mp_size_t rn,
    mp_srcptr m0, mp_srcptr m1, mp_srcptr m2, mp_srcptr m3, mp_size_t mn,
    mp_ptr tp)
{
    if (rn < FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF
        || mn < FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF)
    {
        mp_ptr p0, p1;
        unsigned int i;

        /* Temporary storage: 3 rn + 2 mn */
        p0 = tp + rn;
        p1 = p0 + rn + mn;

        for (i = 0; i < 2; i++)
        {
            flint_mpn_copyi(tp, r0, rn);

            MUL(p0, r0, rn, m0, mn);
            MUL(p1, r1, rn, m3, mn);
            MUL(r0, r1, rn, m2, mn);
            MUL(r1, tp, rn, m1, mn);

            r0[rn + mn] = mpn_add_n(r0, r0, p0, rn + mn);
            r1[rn + mn] = mpn_add_n(r1, r1, p1, rn + mn);

            r0 = r2; r1 = r3;
        }
    }
    else
        matrix22_mul_strassen(r0, r1, r2, r3, rn, m0, m1, m2, m3, mn, tp);
}

/* Multiply M by M1 from the right. Needs 3*(M->n + M1->n) + 5 limbs of
   temporary storage (see _flint_mpn_hgcd_matrix22_mul_strassen_itch). */
void
_flint_mpn_hgcd_matrix_mul(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1, mp_ptr tp)
{
    mp_size_t n;

    /* About the new size of M:s elements. Since M1's diagonal elements are
       > 0, no element can decrease. The new elements are of size
       M->n + M1->n, one limb more or less. The computation of the matrix
       product produces elements of size M->n + M1->n + 1. But the true
       size, after normalization, may be three limbs smaller.

       The reason that the product has normalized size >= M->n + M1->n - 2
       is subtle. It depends on the fact that M and M1 can be factored as
       products of (1,1; 0,1) and (1,0; 1,1), and that we can't have M
       ending with a large power and M1 starting with a large power of the
       same matrix. */

    FLINT_ASSERT(M->n + M1->n < M->alloc);

    FLINT_ASSERT((M->p[0][0][M->n - 1] | M->p[0][1][M->n - 1]
                | M->p[1][0][M->n - 1] | M->p[1][1][M->n - 1]) > 0);

    FLINT_ASSERT((M1->p[0][0][M1->n - 1] | M1->p[0][1][M1->n - 1]
                | M1->p[1][0][M1->n - 1] | M1->p[1][1][M1->n - 1]) > 0);

#if FLINT_HAVE_FFT_SMALL
    if (_flint_mpn_hgcd_matrix_mul_fft(M, M1))
        return;
#endif

    _flint_mpn_hgcd_matrix22_mul_strassen(M->p[0][0], M->p[0][1],
                            M->p[1][0], M->p[1][1], M->n,
                            M1->p[0][0], M1->p[0][1],
                            M1->p[1][0], M1->p[1][1], M1->n, tp);

    /* Index of last potentially non-zero limb, size is one greater. */
    n = M->n + M1->n;

    /* For consecutive half-gcd matrices this loop runs at most three
       times (see above), but other products may shrink more. */
    while (n > 0 && (M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) == 0)
        n--;

    FLINT_ASSERT((M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) > 0);

    M->n = n + 1;
}

/* Multiplies the least significant p limbs of (a;b) by M^-1.
   Temporary space needed: 2 * (p + M->n) */
mp_size_t
_flint_mpn_hgcd_matrix_adjust(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p, mp_ptr tp)
{
    /* M^-1 (a;b) = (r11, -r01; -r10, r00) (a ; b)
       = (r11 a - r01 b; - r10 a + r00 b */

    mp_ptr t0 = tp;
    mp_ptr t1 = tp + p + M->n;
    mp_limb_t ah, bh;
    mp_limb_t cy;

    FLINT_ASSERT(p + M->n < n);

#if FLINT_HAVE_FFT_SMALL
    {
        mp_size_t nn = _flint_mpn_hgcd_matrix_adjust_fft(M, n, ap, bp, p);
        if (nn != 0)
            return nn;
    }
#endif

    /* First compute the two values depending on a, before overwriting a */

    MUL(t0, M->p[1][1], M->n, ap, p);
    MUL(t1, M->p[1][0], M->n, ap, p);

    /* Update a */
    flint_mpn_copyi(ap, t0, p);
    ah = mpn_add(ap + p, ap + p, n - p, t0 + p, M->n);

    MUL(t0, M->p[0][1], M->n, bp, p);

    cy = mpn_sub(ap, ap, n, t0, p + M->n);
    FLINT_ASSERT(cy <= ah);
    ah -= cy;

    /* Update b */
    MUL(t0, M->p[0][0], M->n, bp, p);

    flint_mpn_copyi(bp, t0, p);
    bh = mpn_add(bp + p, bp + p, n - p, t0 + p, M->n);
    cy = mpn_sub(bp, bp, n, t1, p + M->n);
    FLINT_ASSERT(cy <= bh);
    bh -= cy;

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

    Subdivision step and hgcd step

******************************************************************************/

/*
    Used when hgcd or hgcd2 has failed. Then either one of a or b is small,
    or the difference is small. Perform one subtraction followed by one
    division. The normal case is to compute the reduced a and b, and return
    the new size.

    If s == 0 (used for gcd and gcdext), returns zero if the gcd is found.

    If s > 0, don't reduce to size <= s, and return zero if no reduction is
    possible (if either a, b or |a-b| is of size <= s).

    The hook function is called as hook(ctx, gp, gn, qp, qn, d) in the
    following cases:

    + If A = B at the start, G is the gcd, Q is NULL, d = -1.

    + If one input is zero at the start, G is the gcd, Q is NULL, d = 0 if
      A = G and d = 1 if B = G.

    Otherwise, if d = 0 we have just subtracted a multiple of A from B, and
    if d = 1 we have subtracted a multiple of B from A.

    + If A = B after subtraction, G is the gcd, Q is NULL.

    + If we get a zero remainder after division, G is the gcd, Q is the
      quotient.

    + Otherwise, G is NULL, Q is the quotient (often 1).
*/
mp_size_t
_flint_mpn_gcd_subdiv_step(mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t s,
    _flint_mpn_gcd_subdiv_step_hook * hook, void * ctx, mp_ptr tp)
{
    static const mp_limb_t one = 1;
    mp_size_t an, bn, qn;
    int swapped;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(ap[n - 1] > 0 || bp[n - 1] > 0);

    an = bn = n;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);

    swapped = 0;

    /* Arrange so that a < b, subtract b -= a, and maintain normalization. */
    if (an == bn)
    {
        int c = mpn_cmp(ap, bp, an);
        if (FLINT_UNLIKELY(c == 0))
        {
            /* For gcdext, return the smallest of the two cofactors, so
               pass d = -1. */
            if (s == 0)
                hook(ctx, ap, an, NULL, 0, -1);
            return 0;
        }
        else if (c > 0)
        {
            FLINT_SWAP(mp_ptr, ap, bp);
            swapped ^= 1;
        }
    }
    else
    {
        if (an > bn)
        {
            MPN_SWAP(ap, an, bp, bn);
            swapped ^= 1;
        }
    }
    if (an <= s)
    {
        if (s == 0)
            hook(ctx, bp, bn, NULL, 0, swapped ^ 1);
        return 0;
    }

    mpn_sub(bp, bp, bn, ap, an);
    MPN_NORM(bp, bn);
    FLINT_ASSERT(bn > 0);

    if (bn <= s)
    {
        /* Undo subtraction. */
        mp_limb_t cy = mpn_add(bp, ap, an, bp, bn);
        if (cy > 0)
            bp[an] = cy;
        return 0;
    }

    /* Arrange so that a < b */
    if (an == bn)
    {
        int c = mpn_cmp(ap, bp, an);
        if (FLINT_UNLIKELY(c == 0))
        {
            if (s > 0)
                /* Just record subtraction and return */
                hook(ctx, NULL, 0, &one, 1, swapped);
            else
                /* Found gcd. */
                hook(ctx, bp, bn, NULL, 0, swapped);
            return 0;
        }

        hook(ctx, NULL, 0, &one, 1, swapped);

        if (c > 0)
        {
            FLINT_SWAP(mp_ptr, ap, bp);
            swapped ^= 1;
        }
    }
    else
    {
        hook(ctx, NULL, 0, &one, 1, swapped);

        if (an > bn)
        {
            MPN_SWAP(ap, an, bp, bn);
            swapped ^= 1;
        }
    }

    flint_mpn_tdiv_qr(tp, bp, bp, bn, ap, an);
    qn = bn - an + 1;
    bn = an;
    MPN_NORM(bp, bn);

    if (FLINT_UNLIKELY(bn <= s))
    {
        if (s == 0)
        {
            hook(ctx, ap, an, tp, qn, swapped);
            return 0;
        }

        /* Quotient is one too large, so decrement it and add back A. */
        if (bn > 0)
        {
            mp_limb_t cy = mpn_add(bp, ap, an, bp, bn);
            if (cy)
                bp[an++] = cy;
        }
        else
            flint_mpn_copyi(bp, ap, an);

        mpn_sub_1(tp, tp, qn, 1);
    }
    hook(ctx, NULL, 0, tp, qn, swapped);
    return an;
}

static void
hgcd_hook(void * p, mp_srcptr gp, mp_size_t gn, mp_srcptr qp, mp_size_t qn, int d)
{
    FLINT_ASSERT(!gp);
    FLINT_ASSERT(d >= 0);
    FLINT_ASSERT(d <= 1);

    (void) gp;
    (void) gn;

    MPN_NORM(qp, qn);
    if (qn > 0)
    {
        flint_mpn_hgcd_matrix_struct * M = (flint_mpn_hgcd_matrix_struct *) p;
        /* NOTES: This is a bit ugly. A tp area is passed to
           gcd_subdiv_step, which stores q at the start of that area. We
           now use the rest. */
        mp_ptr tp = (mp_ptr) qp + qn;
        _flint_mpn_hgcd_matrix_update_q(M, qp, qn, d, tp);
    }
}

/*
    Perform a few steps, using some of hgcd2, subtraction and division.
    Reduces the size by almost one limb or more, but never below the given
    size s. Return new size for a and b, or 0 if no more steps are possible.

    If hgcd2 succeeds, needs temporary space for hgcd_matrix_mul_1, M->n
    limbs, and hgcd_mul_matrix1_inverse_vector, n limbs. If hgcd2 fails,
    needs space for the quotient, qn <= n - s limbs, for and
    hgcd_matrix_update_q, qn + (size of the appropriate column of M) <=
    (resulting size of M) + 1.

    If N is the input size to the calling hgcd, then s = floor(N/2) + 1,
    M->n < N, qn + product size <= n - s + n - s + 1 = 2 (n - s) + 1 <= N.
*/
mp_size_t
_flint_mpn_hgcd_step(mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t s,
    flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
{
    flint_mpn_hgcd_matrix1_struct M1;
    mp_limb_t mask;
    mp_limb_t ah, al, bh, bl;

    FLINT_ASSERT(n > s);

    mask = ap[n - 1] | bp[n - 1];
    FLINT_ASSERT(mask > 0);

    if (n == s + 1)
    {
        if (mask < 4)
            goto subtract;

        ah = ap[n - 1]; al = ap[n - 2];
        bh = bp[n - 1]; bl = bp[n - 2];
    }
    else if (mask & (UWORD(1) << (FLINT_BITS - 1)))
    {
        ah = ap[n - 1]; al = ap[n - 2];
        bh = bp[n - 1]; bl = bp[n - 2];
    }
    else
    {
        int shift = flint_clz(mask);
        ah = EXTRACT_LIMB(shift, ap[n - 1], ap[n - 2]);
        al = EXTRACT_LIMB(shift, ap[n - 2], ap[n - 3]);
        bh = EXTRACT_LIMB(shift, bp[n - 1], bp[n - 2]);
        bl = EXTRACT_LIMB(shift, bp[n - 2], bp[n - 3]);
    }

    /* Try an hgcd2 step */
    if (_flint_mpn_hgcd2(ah, al, bh, bl, &M1))
    {
        /* Multiply M <- M * M1 */
        _flint_mpn_hgcd_matrix_mul_1(M, &M1, tp);

        /* Multiply M1^{-1} (a;b), in place */
        return _flint_mpn_hgcd_mul_matrix1_inverse_vector(&M1, ap, ap, bp, n);
    }

subtract:
    return _flint_mpn_gcd_subdiv_step(ap, bp, n, s, hgcd_hook, M, tp);
}

/******************************************************************************

    hgcd

******************************************************************************/

/*
    Size analysis for hgcd:

    For the recursive calls, we have n1 <= ceil(n / 2). Then the storage need
    is determined by the storage for the recursive call computing M1, and
    hgcd_matrix_adjust and hgcd_matrix_mul calls that use M1 (after this, the
    storage needed for M1 can be recycled).

    Let S(r) denote the required storage. For M1 we need 4 * (ceil(n1/2) + 1)
    = 4 * (ceil(n/4) + 1), for the hgcd_matrix_adjust call, we need n + 2,
    and for the hgcd_matrix_mul, we may need 3 ceil(n/2) + 8. In total,
    4 * ceil(n/4) + 3 ceil(n/2) + 12 <= 10 ceil(n/4) + 12.

    For the recursive call, we need S(n1) = S(ceil(n/2)).

    S(n) <= 10*ceil(n/4) + 12 + S(ceil(n/2))
         <= 10*(ceil(n/4) + ... + ceil(n/2^(1+k))) + 12k + S(ceil(n/2^k))
         <= 10*(2 ceil(n/4) + k) + 12k + S(ceil(n/2^k))
         <= 20 ceil(n/4) + 22k + S(ceil(n/2^k))
*/

mp_size_t
_flint_mpn_hgcd_itch(mp_size_t n)
{
    unsigned int k;
    mp_size_t nscaled;

    if (n < FLINT_MPN_HGCD_CUTOFF)
        return n;

    /* Get the recursion depth. */
    nscaled = (n - 1) / (FLINT_MPN_HGCD_CUTOFF - 1);
    k = FLINT_BITS - flint_clz(nscaled);

    return 20 * ((n + 3) / 4) + 22 * k + FLINT_MPN_HGCD_CUTOFF;
}

/* Reduces a,b until |a-b| fits in n/2 + 1 limbs. Constructs matrix M with
   elements of size at most (n+1)/2 - 1. Returns new size of a, b, or zero
   if no reduction is possible. */
mp_size_t
_flint_mpn_hgcd(mp_ptr ap, mp_ptr bp, mp_size_t n,
    flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
{
    mp_size_t s = n / 2 + 1;
    mp_size_t nn;
    int success = 0;

    if (n <= s)
        /* Happens when n <= 2, a fairly uninteresting case but exercised by
           the random inputs of the testsuite. */
        return 0;

    FLINT_ASSERT((ap[n - 1] | bp[n - 1]) > 0);

    FLINT_ASSERT((n + 1) / 2 - 1 < M->alloc);

    if (n >= FLINT_MPN_HGCD_CUTOFF)
    {
        mp_size_t n2 = (3 * n) / 4 + 1;
        mp_size_t p = n / 2;

        nn = _flint_mpn_hgcd_reduce(M, ap, bp, n, p, tp);
        if (nn)
        {
            n = nn;
            success = 1;
        }

        /* NOTE: It appears this loop never runs more than once (at least
           when not recursing to hgcd_appr). */
        while (n > n2)
        {
            /* Needs n + 1 storage */
            nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);
            if (!nn)
                return success ? n : 0;

            n = nn;
            success = 1;
        }

        if (n > s + 2)
        {
            flint_mpn_hgcd_matrix_struct M1;
            mp_size_t scratch;

            p = 2 * s - n + 1;
            scratch = FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n - p);

            _flint_mpn_hgcd_matrix_init(&M1, n - p, tp);

            nn = _flint_mpn_hgcd(ap + p, bp + p, n - p, &M1, tp + scratch);
            if (nn > 0)
            {
                /* We always have max(M) > 2^{-(FLINT_BITS + 1)} max(M1) */
                FLINT_ASSERT(M->n + 2 >= M1.n);

                /* Furthermore, assume M ends with a quotient (1, q; 0, 1),
                   then either q or q + 1 is a correct quotient, and M1 will
                   start with either (1, 0; 1, 1) or (2, 1; 1, 1). This rules
                   out the case that the size of M * M1 is much smaller than
                   the expected M->n + M1->n. */

                FLINT_ASSERT(M->n + M1.n < M->alloc);

#if FLINT_HAVE_FFT_SMALL
                /* both operations with one transform of M1 */
                {
                    mp_size_t n2;
                    n2 = _flint_mpn_hgcd_matrix_adjust_mul_fft(M, &M1, p + nn, ap, bp, p);
                    if (n2 != 0)
                    {
                        n = n2;
                        success = 1;
                        goto second_done;
                    }
                }
#endif

                /* Needs 2 (p + M->n) <= 2 (2*s - n2 + 1 + n2 - s - 1)
                   = 2*s <= 2*(floor(n/2) + 1) <= n + 2. */
                n = _flint_mpn_hgcd_matrix_adjust(&M1, p + nn, ap, bp, p, tp + scratch);

                /* We need a bound for of M->n + M1.n. Let n be the original
                   input size. Then

                   ceil(n/2) - 1 >= size of product >= M.n + M1.n - 2

                   and it follows that

                   M.n + M1.n <= ceil(n/2) + 1

                   Then 3*(M.n + M1.n) + 5 <= 3 * ceil(n/2) + 8 is the amount
                   of needed scratch space. */
                _flint_mpn_hgcd_matrix_mul(M, &M1, tp + scratch);
                success = 1;
            }
        }
    }

#if FLINT_HAVE_FFT_SMALL
second_done:
#endif
    for (;;)
    {
        /* Needs s+3 < n */
        nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);
        if (!nn)
            return success ? n : 0;

        n = nn;
        success = 1;
    }
}

/******************************************************************************

    hgcd_appr

******************************************************************************/

/* Identical to _flint_mpn_hgcd_itch. */
mp_size_t
_flint_mpn_hgcd_appr_itch(mp_size_t n)
{
    if (n < FLINT_MPN_HGCD_APPR_CUTOFF)
        return n;
    else
    {
        unsigned int k;
        mp_size_t nscaled;

        /* Get the recursion depth. */
        nscaled = (n - 1) / (FLINT_MPN_HGCD_APPR_CUTOFF - 1);
        k = FLINT_BITS - flint_clz(nscaled);

        return 20 * ((n + 3) / 4) + 22 * k + FLINT_MPN_HGCD_CUTOFF;
    }
}

/* Destroys inputs. */
int
_flint_mpn_hgcd_appr(mp_ptr ap, mp_ptr bp, mp_size_t n,
    flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
{
    mp_size_t s;
    int success = 0;

    FLINT_ASSERT(n > 0);

    FLINT_ASSERT((ap[n - 1] | bp[n - 1]) != 0);

    if (n <= 2)
        /* Implies s = n. A fairly uninteresting case but exercised by the
           random inputs of the testsuite. */
        return 0;

    FLINT_ASSERT((n + 1) / 2 - 1 < M->alloc);

    /* We aim for reduction of to FLINT_BITS * s bits. But each time we
       discard some of the least significant limbs, we must keep one
       additional bit to account for the truncation error. We maintain the
       FLINT_BITS * s - extra_bits as the current target size. */

    s = n / 2 + 1;
    if (n < FLINT_MPN_HGCD_APPR_CUTOFF)
    {
        unsigned int extra_bits = 0;

        while (n > 2)
        {
            mp_size_t nn;

            FLINT_ASSERT(n > s);
            FLINT_ASSERT(n <= 2 * s);

            nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);
            if (!nn)
                break;

            n = nn;
            success = 1;

            /* We can truncate and discard the lower p bits whenever
               nbits <= 2*sbits - p. To account for the truncation error, we
               must adjust

               sbits <-- sbits + 1 - p,

               rather than just sbits <-- sbits - p. This adjustment makes
               the produced matrix slightly smaller than it could be. */

            if (FLINT_BITS * (n + 1) + 2 * extra_bits <= 2 * FLINT_BITS * s)
            {
                mp_size_t p = (FLINT_BITS * (2 * s - n) - 2 * extra_bits) / FLINT_BITS;

                if (extra_bits == 0)
                {
                    /* We cross a limb boundary and bump s. We can't do that
                       if the result is that it makes makes min(U, V)
                       smaller than 2^{FLINT_BITS} s. */
                    if (s + 1 == n
                        || flint_mpn_zero_p(ap + s + 1, n - s - 1)
                        || flint_mpn_zero_p(bp + s + 1, n - s - 1))
                        continue;

                    extra_bits = FLINT_BITS - 1;
                    s++;
                }
                else
                {
                    extra_bits--;
                }

                /* Drop the p least significant limbs */
                ap += p; bp += p; n -= p; s -= p;
            }
        }

        FLINT_ASSERT(s > 0);

        if (extra_bits > 0)
        {
            /* We can get here only of we have dropped at least one of the
               least significant bits, so we can decrement ap and bp. We can
               then shift left extra bits using mpn_rshift. */
            /* NOTE: In the unlikely case that n is large, it would be
               preferable to do an initial subdiv step to reduce the size
               before shifting, but that would mean duplicating
               gcd_subdiv_step with a bit count rather than a limb count. */
            ap--; bp--;
            ap[0] = mpn_rshift(ap + 1, ap + 1, n, FLINT_BITS - extra_bits);
            bp[0] = mpn_rshift(bp + 1, bp + 1, n, FLINT_BITS - extra_bits);
            n += (ap[n] | bp[n]) > 0;

            FLINT_ASSERT(success);

            while (n > 2)
            {
                mp_size_t nn;

                FLINT_ASSERT(n > s);
                FLINT_ASSERT(n <= 2 * s);

                nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);

                if (!nn)
                    return 1;

                n = nn;
            }
        }

        if (n == 2)
        {
            flint_mpn_hgcd_matrix1_struct M1;
            FLINT_ASSERT(s == 1);

            if (_flint_mpn_hgcd2(ap[1], ap[0], bp[1], bp[0], &M1))
            {
                /* Multiply M <- M * M1 */
                _flint_mpn_hgcd_matrix_mul_1(M, &M1, tp);
                success = 1;
            }
        }
        return success;
    }
    else
    {
        mp_size_t n2 = (3 * n) / 4 + 1;
        mp_size_t p = n / 2;
        mp_size_t nn;

        nn = _flint_mpn_hgcd_reduce(M, ap, bp, n, p, tp);
        if (nn)
        {
            n = nn;
            /* FIXME: Discard some of the low limbs immediately? */
            success = 1;
        }

        while (n > n2)
        {
            /* Needs n + 1 storage */
            nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);
            if (!nn)
                return success;

            n = nn;
            success = 1;
        }
        if (n > s + 2)
        {
            flint_mpn_hgcd_matrix_struct M1;
            mp_size_t scratch;

            p = 2 * s - n + 1;
            scratch = FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n - p);

            _flint_mpn_hgcd_matrix_init(&M1, n - p, tp);
            if (_flint_mpn_hgcd_appr(ap + p, bp + p, n - p, &M1, tp + scratch))
            {
                /* We always have max(M) > 2^{-(FLINT_BITS + 1)} max(M1) */
                FLINT_ASSERT(M->n + 2 >= M1.n);

                /* Furthermore, assume M ends with a quotient (1, q; 0, 1),
                   then either q or q + 1 is a correct quotient, and M1 will
                   start with either (1, 0; 1, 1) or (2, 1; 1, 1). This rules
                   out the case that the size of M * M1 is much smaller than
                   the expected M->n + M1->n. */

                FLINT_ASSERT(M->n + M1.n < M->alloc);

                /* We need a bound for of M->n + M1.n. Let n be the original
                   input size. Then

                   ceil(n/2) - 1 >= size of product >= M.n + M1.n - 2

                   and it follows that

                   M.n + M1.n <= ceil(n/2) + 1

                   Then 3*(M.n + M1.n) + 5 <= 3 * ceil(n/2) + 8 is the amount
                   of needed scratch space. */
                _flint_mpn_hgcd_matrix_mul(M, &M1, tp + scratch);
                return 1;
            }
        }

        for (;;)
        {
            mp_size_t nn;

            FLINT_ASSERT(n > s);
            FLINT_ASSERT(n <= 2 * s);

            nn = _flint_mpn_hgcd_step(n, ap, bp, s, M, tp);

            if (!nn)
                return success;

            n = nn;
            success = 1;
        }
    }
}

/******************************************************************************

    hgcd_reduce

******************************************************************************/

/* Computes R -= A * B. Result must be non-negative. Normalized down to size
   an, and resulting size is returned. */
static mp_size_t
submul(mp_ptr rp, mp_size_t rn, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_ptr tp;
    TMP_INIT;

    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(rn >= an);
    FLINT_ASSERT(an + bn <= rn + 1);

    TMP_START;
    tp = TMP_ALLOC((an + bn) * sizeof(mp_limb_t));

    flint_mpn_mul(tp, ap, an, bp, bn);
    FLINT_ASSERT((an + bn <= rn) || (tp[rn] == 0));
    mpn_sub(rp, rp, rn, tp, an + bn - (an + bn > rn));
    TMP_END;

    while (rn > an && (rp[rn - 1] == 0))
        rn--;

    return rn;
}

/*
    Computes (a, b) <-- M^{-1} (a; b), using that the result is known to be
    small: both entries fit in nn = max(un, vn) limbs (see below), so it
    suffices to compute them modulo B^nn (or modulo B^N - 1 for N > nn as in
    GMP, which uses mpn_mulmod_bnm1). The fft_small version computes the
    four products in transform space; otherwise low products are used.
*/
static mp_size_t
hgcd_matrix_apply(const flint_mpn_hgcd_matrix_struct * M, mp_ptr ap, mp_ptr bp, mp_size_t n)
{
    mp_size_t an, bn, un, vn, nn;
    mp_size_t mn[2][2];
    unsigned int i, j;

    FLINT_ASSERT((ap[n - 1] | bp[n - 1]) > 0);

    an = n;
    MPN_NORM(ap, an);
    bn = n;
    MPN_NORM(bp, bn);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            mp_size_t k;
            k = M->n;
            MPN_NORM(M->p[i][j], k);
            mn[i][j] = k;
        }

    FLINT_ASSERT(mn[0][0] > 0);
    FLINT_ASSERT(mn[1][1] > 0);
    FLINT_ASSERT((mn[0][1] | mn[1][0]) > 0);

    if (mn[0][1] == 0)
    {
        /* A unchanged, M = (1, 0; q, 1) */
        FLINT_ASSERT(mn[0][0] == 1);
        FLINT_ASSERT(M->p[0][0][0] == 1);
        FLINT_ASSERT(mn[1][1] == 1);
        FLINT_ASSERT(M->p[1][1][0] == 1);

        /* Put B <-- B - q A */
        nn = submul(bp, bn, ap, an, M->p[1][0], mn[1][0]);
    }
    else if (mn[1][0] == 0)
    {
        /* B unchanged, M = (1, q; 0, 1) */
        FLINT_ASSERT(mn[0][0] == 1);
        FLINT_ASSERT(M->p[0][0][0] == 1);
        FLINT_ASSERT(mn[1][1] == 1);
        FLINT_ASSERT(M->p[1][1][0] == 1);

        /* Put A  <-- A - q * B */
        nn = submul(ap, an, bp, bn, M->p[0][1], mn[0][1]);
    }
    else
    {
        /* A = m00 a + m01 b  ==> a <= A / m00, b <= A / m01.
           B = m10 a + m11 b  ==> a <= B / m10, b <= B / m11. */
        un = FLINT_MIN(an - mn[0][0], bn - mn[1][0]) + 1;
        vn = FLINT_MIN(an - mn[0][1], bn - mn[1][1]) + 1;

        nn = FLINT_MAX(un, vn);

#if FLINT_HAVE_FFT_SMALL
        if (!_flint_mpn_hgcd_matrix_apply_fft(M, mn, ap, bp, n, nn))
#endif
        {
            mp_ptr tp, sp;
            mp_size_t k;
            TMP_INIT;

            TMP_START;
            tp = TMP_ALLOC(2 * (nn + M->n) * sizeof(mp_limb_t));
            sp = tp + nn + M->n;

            /* Only the low nn limbs are needed. */
            k = FLINT_MIN(n, nn);

            /* new a = m11 a - m01 b, new b = m00 b - m10 a (mod B^nn) */
            flint_mpn_zero(tp, nn + M->n);
            flint_mpn_zero(sp, nn + M->n);
            MUL(tp, ap, k, M->p[1][1], FLINT_MIN(mn[1][1], nn));
            MUL(sp, bp, k, M->p[0][1], FLINT_MIN(mn[0][1], nn));
            mpn_sub_n(tp, tp, sp, nn);

            MUL(sp, ap, k, M->p[1][0], FLINT_MIN(mn[1][0], nn));
            flint_mpn_copyi(ap, tp, nn);
            flint_mpn_zero(tp, nn + M->n);
            MUL(tp, bp, k, M->p[0][0], FLINT_MIN(mn[0][0], nn));
            mpn_sub_n(tp, tp, sp, nn);
            flint_mpn_copyi(bp, tp, nn);

            TMP_END;
        }

        while ((ap[nn - 1] | bp[nn - 1]) == 0)
        {
            nn--;
            FLINT_ASSERT(nn > 0);
        }
    }

    return nn;
}

mp_size_t
_flint_mpn_hgcd_reduce_itch(mp_size_t n, mp_size_t p)
{
    mp_size_t itch;
    if (n < FLINT_MPN_HGCD_REDUCE_CUTOFF)
    {
        itch = _flint_mpn_hgcd_itch(n - p);

        /* For arbitrary p, the storage for _adjust is 2*(p + M->n) = 2 *
           (p + ceil((n-p)/2) - 1 <= n + p - 1 */
        if (itch < n + p - 1)
            itch = n + p - 1;
    }
    else
    {
        itch = 2 * (n - p) + _flint_mpn_hgcd_itch(n - p);
        /* Currently, hgcd_matrix_apply allocates its own storage. */
    }
    return itch;
}

mp_size_t
_flint_mpn_hgcd_reduce(flint_mpn_hgcd_matrix_struct * M,
    mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t p, mp_ptr tp)
{
    mp_size_t nn;
    if (n < FLINT_MPN_HGCD_REDUCE_CUTOFF)
    {
        nn = _flint_mpn_hgcd(ap + p, bp + p, n - p, M, tp);
        if (nn > 0)
            /* Needs 2*(p + M->n) <= 2*(floor(n/2) + ceil(n/2) - 1)
               = 2 (n - 1) */
            return _flint_mpn_hgcd_matrix_adjust(M, p + nn, ap, bp, p, tp);
    }
    else
    {
        flint_mpn_copyi(tp, ap + p, n - p);
        flint_mpn_copyi(tp + n - p, bp + p, n - p);
        if (_flint_mpn_hgcd_appr(tp, tp + n - p, n - p, M, tp + 2 * (n - p)))
            return hgcd_matrix_apply(M, ap, bp, n);
    }
    return 0;
}
