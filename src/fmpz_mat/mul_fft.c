/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fft.h"

/*
    fft coeffs are mod m = 2^(FLINT_BITS*limbs) + 1.
    Viewing these coefficients as in the signed range (-m/2, m/2):
        set the mpn z to the twos complement evaluation at 2^bits.
    The lower zn limbs are written to z and the lowest bit of the
    nominal zn^th limb in the twos complement representation is returned.

    The behaviour of this function does NOT depend on the initial value of z.
*/
static mp_limb_t fft_combine_bits_signed(
    mp_limb_t * z,
    mp_limb_t ** a, mp_size_t alen,
    flint_bitcnt_t bits,
    mp_size_t limbs,
    mp_size_t zn)
{
    mp_size_t i, zout;
    mp_limb_t * t;
    mp_limb_t f;
    TMP_INIT;

    FLINT_ASSERT(bits > 1);

    TMP_START;

    t = TMP_ARRAY_ALLOC((limbs + 1), mp_limb_t);

    f = 0;
    zout = 0;

    for (i = 0; i < alen; i++)
    {
        /* add the i^th coeffs a[i] */
        mp_limb_t q = (bits*i)/FLINT_BITS;
        mp_limb_t r = (bits*i)%FLINT_BITS;
        mp_limb_t s;
        mp_limb_t halflimb = UWORD(1) << (FLINT_BITS - 1);

        if (a[i][limbs] | (a[i][limbs - 1] > halflimb))
        {
            mpn_sub_1(t, a[i], limbs, UWORD(1));
            s = 1;
        }
        else
        {
            mpn_copyi(t, a[i], limbs);
            s = 0;
        }

        t[limbs] = -s;

        /*
            t[0] ... t[limbs] is now the twos complement signed version of a[i]
            and the shift occurs without overflow
        */

        if (r != 0)
            mpn_lshift(t, t, limbs + 1, r);

        if (q < zn)
        {
            size_t new_zout = FLINT_MIN(zn, q + limbs + 1);
            FLINT_ASSERT(new_zout >= zout);

            while (zout < new_zout)
                z[zout++] = -f;

            FLINT_ASSERT(new_zout > q);
            f ^= s;
            f ^= mpn_add_n(z + q, z + q, t, new_zout - q);
        }
        else
        {
            if (q == zn)
                f ^= t[0]&1;
            break;
        }
    }

    while (zout < zn)
        z[zout++] = -f;

    TMP_END;

    FLINT_ASSERT(f == 0 || f == 1);
    return f;
}

/*
    Split into coefficients from |x| evaluated at 2^bits,
    and do a negmod on each coefficient for x < 0.
*/
static mp_size_t fft_split_bits_fmpz(
    mp_limb_t ** poly,
    const fmpz_t x,
    flint_bitcnt_t bits,
    mp_size_t limbs)
{
    mp_size_t len;
    int x_is_neg = 0;

    if (COEFF_IS_MPZ(*x))
    {
        mp_size_t s = COEFF_TO_PTR(*x)->_mp_size;
        x_is_neg = s < 0;
        len = fft_split_bits(poly, COEFF_TO_PTR(*x)->_mp_d,
                             x_is_neg ? -s : s, bits, limbs);
    }
    else if (!fmpz_is_zero(x))
    {
        mp_limb_t ux;
        x_is_neg = *x < 0;
        ux = x_is_neg ? -*x : *x;
        len = fft_split_bits(poly, &ux, 1, bits, limbs);
    }
    else
    {
        len = 0;
    }

    if (x_is_neg)
    {
        mp_size_t i;
        for (i = 0; i < len; i++)
            mpn_negmod_2expp1(poly[i], poly[i], limbs);
    }

    return len;
}

static void fft_combine_bits_fmpz(
    fmpz_t x,
    mp_limb_t ** poly, slong length,
    flint_bitcnt_t bits,
    mp_size_t limbs,
    mp_size_t total_limbs,
    int sign)
{
    __mpz_struct * mx = _fmpz_promote(x);
    mp_limb_t * d = FLINT_MPZ_REALLOC(mx, total_limbs);
    if (sign)
    {
        if (fft_combine_bits_signed(d, poly, length, bits, limbs, total_limbs))
        {
            mpn_neg(d, d, total_limbs);
            MPN_NORM(d, total_limbs);
            /* total_limbs should have started high enough to prevent d = 0 here */
            FLINT_ASSERT(total_limbs > 0);
            mx->_mp_size = -total_limbs;
        }
        else
        {
            MPN_NORM(d, total_limbs);
            mx->_mp_size = total_limbs;
        }
    }
    else
    {
        flint_mpn_zero(d, total_limbs);
        fft_combine_bits(d, poly, length, bits, limbs, total_limbs);
        MPN_NORM(d, total_limbs);
        mx->_mp_size = total_limbs;
    }
    _fmpz_demote_val(x);
}


static void _either_fft_or_mfa(
    ulong ** coeffs,
    slong n, flint_bitcnt_t w,
    ulong ** t1, ulong ** t2, ulong ** t3,
    slong n1,
    flint_bitcnt_t depth,
    slong trunc,
    slong limbs,
    int use_mfa)
{
    ulong trunc2, rs, s, u;
    slong l;

    if (use_mfa)
    {
        fft_mfa_truncate_sqrt2(coeffs, n, w, t1, t2, t3, n1, trunc);

        for (l = 0; l < 2*n; l++)
            mpn_normmod_2expp1(coeffs[l], limbs);

        /* the second half is out of order */
        trunc2 = (trunc - 2*n)/n1;
        for (s = 0; s < trunc2; s++)
        {
            rs = n_revbin(s, depth - depth/2 + 1);
            for (u = 0; u < n1; u++)
            {
                l = 2*n + rs*n1 + u;
                mpn_normmod_2expp1(coeffs[l], limbs);
            }
        }
    }
    else
    {
        fft_truncate_sqrt2(coeffs, n, w, t1, t2, t3, trunc);

        for (l = 0; l < trunc; l++)
            mpn_normmod_2expp1(coeffs[l], limbs);
    }
}

static void _dot(
    ulong* c,
    ulong** A, slong Astride,
    ulong** B, slong Bstride,
    slong len,
    slong limbs,
    ulong* t,   /* temp of length limbs + 1 */
    ulong* t2)  /* temp of length limbs + 1 */
{
    slong i;
    flint_bitcnt_t nw = limbs*FLINT_BITS; /* n*w */

    FLINT_ASSERT(len > 0);

    i = 0;
    do {
        const ulong* a = A[i*Astride];
        const ulong* b = B[i*Bstride];
        if (i == 0)
        {
            c[limbs] = flint_mpn_mulmod_2expp1_basecase(c, a, b,
                                                 2*a[limbs] + b[limbs], nw, t2);
        }
        else
        {
            t[limbs] = flint_mpn_mulmod_2expp1_basecase(t, a, b,
                                                 2*a[limbs] + b[limbs], nw, t2);

            c[limbs] += t[limbs];
            c[limbs] += mpn_add_n(c, c, t, limbs);

            mpn_normmod_2expp1(c, limbs);
        }
    } while (++i < len);
}


void _fmpz_mat_mul_truncate_sqrt2(
    fmpz_mat_t C,
    const fmpz_mat_t A, slong Abits,
    const fmpz_mat_t B, slong Bbits,
    flint_bitcnt_t depth,
    flint_bitcnt_t w,
    slong j1, slong j2,
    const int use_mfa,
    const int sign)
{
    slong M = fmpz_mat_nrows(A);
    slong K = fmpz_mat_ncols(A);
    slong N = fmpz_mat_ncols(B);
    slong clgK = FLINT_CLOG2(K) + sign;
    slong n = WORD(1) << depth;
    ulong trunc, sqrt;
    ulong bits1 = (n*w - (depth + 1 + clgK))/2;
/* This bound on Climbs is also ok with the original clgK = FLINT_CLOG2(K) */
    slong Climbs = (Abits + Bbits + clgK + FLINT_BITS - 1)/FLINT_BITS;
    slong limbs = (n*w)/FLINT_BITS;
    slong size = limbs + 1;
    slong i, j, l, h;
    ulong * temp, *t, * t1, * t2, * t3, * Adata, * Bdata, * Cdata;
    ulong ** coeffs, ** Acoeffs, ** Bcoeffs, ** Ccoeffs;

    FLINT_ASSERT(limbs > 0);
    FLINT_ASSERT(limbs*FLINT_BITS == n*w);

    /* should have K*min(j1,j2)*2^bits1 <= 2^(n*w) */
    FLINT_ASSERT(j1 <= 2*n || j2 <= 2*n);

    /*
        1 array of length 2*size for mulmod_2expp1.
        1 array of length size for addmulmod
        3 arrays of length size for fft/ifft.
        M*K*4*n arrays of length size for A's fft rep
        K*N*4*n arrays of length size for B's fft rep
            4*n arrays of length size for C's fft rep
    */
    temp = FLINT_ARRAY_ALLOC((6 + 4*n*(M*K + K*N + 1))*size, mp_limb_t);
    t = temp + 2*size;
    t1 = t + size;
    t2 = t1 + size;
    t3 = t2 + size;
    Adata = t3 + size;
    Bdata = Adata + 4*n*M*K*size;
    Cdata = Bdata + 4*n*K*N*size;

    /*
        M*K arrays of pointers of length 4*n for A's coeffs
        K*N arrays of pointers of length 4*n for B's coeffs
          1 array  of pointers of length 4*n for C's coeffs
    */
    coeffs = FLINT_ARRAY_ALLOC(4*n*(M*K + K*N + 1), mp_limb_t*);
    Acoeffs = coeffs;
    Bcoeffs = Acoeffs + 4*n*M*K;
    Ccoeffs = Bcoeffs + 4*n*K*N;

    /* fill in coefficient pointers */
    for (i = 0; i < M; i++)
    for (j = 0; j < K; j++)
    for (l = 0; l < 4*n; l++)
        Acoeffs[(i*K + j)*4*n + l] = Adata + ((i*K + j)*4*n + l)*size;

    for (i = 0; i < K; i++)
    for (j = 0; j < N; j++)
    for (l = 0; l < 4*n; l++)
        Bcoeffs[(i*N + j)*4*n + l] = Bdata + ((i*N + j)*4*n + l)*size;

    for (l = 0; l < 4*n; l++)
        Ccoeffs[l] = Cdata + l*size;

    /* chop inputs into bits */
    for (i = 0; i < M; i++)
    for (j = 0; j < K; j++)
    {
        h = fft_split_bits_fmpz(Acoeffs + (i*K + j)*4*n,
                                fmpz_mat_entry(A, i, j), bits1, limbs);
        for (l = h; l < 4*n; l++)
            flint_mpn_zero(Acoeffs[(i*K + j)*4*n + l], size);
    }

    for (i = 0; i < K; i++)
    for (j = 0; j < N; j++)
    {
        h = fft_split_bits_fmpz(Bcoeffs + (i*N + j)*4*n,
                                fmpz_mat_entry(B, i, j), bits1, limbs);
        for (l = h; l < 4*n; l++)
            flint_mpn_zero(Bcoeffs[(i*N + j)*4*n + l], size);
    }

    FLINT_ASSERT(j1 > 0);
    FLINT_ASSERT(j2 > 0);
    FLINT_ASSERT(j1 + j2 - 1 <= 4*n);

    trunc = j1 + j2 - 1;
    trunc = FLINT_MAX(trunc, 2*n + 1);
    if (use_mfa)
    {
        sqrt = UWORD(1) << (depth/2);
        trunc = (trunc + 2*sqrt - 1) & (-2*sqrt);
    }
    else
    {
        sqrt = 1;
        trunc = (trunc + 1) & -UWORD(2);
    }

    FLINT_ASSERT(trunc > 2*n);
    FLINT_ASSERT(trunc % (2*sqrt) == 0);

    for (i = 0; i < M; i++)
    for (j = 0; j < K; j++)
    {
        _either_fft_or_mfa(Acoeffs + (i*K + j)*4*n, n, w,
                            &t1, &t2, &t3, sqrt, depth, trunc, limbs, use_mfa);
    }

    for (i = 0; i < K; i++)
    for (j = 0; j < N; j++)
    {
        _either_fft_or_mfa(Bcoeffs + (i*N + j)*4*n, n, w,
                            &t1, &t2, &t3, sqrt, depth, trunc, limbs, use_mfa);
    }

    /* pointwise dot products for C[i,j] */
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++)
    {
        if (use_mfa)
        {
            ulong trunc2, rs, s, u;

            for (l = 0; l < 2*n; l++)
            {
                _dot(Ccoeffs[l], Acoeffs + (i*K + 0)*4*n + l, 4*n,
                                 Bcoeffs + (0*N + j)*4*n + l, N*4*n,
                     K, limbs, t, temp);
            }

            /* second half is out of order */
            trunc2 = (trunc - 2*n)/sqrt;
            for (s = 0; s < trunc2; s++)
            {
                rs = n_revbin(s, depth - depth/2 + 1);
                for (u = 0; u < sqrt; u++)
                {
                    l = 2*n + rs*sqrt + u;

                    _dot(Ccoeffs[l], Acoeffs + (i*K + 0)*4*n + l, 4*n,
                                     Bcoeffs + (0*N + j)*4*n + l, N*4*n,
                         K, limbs, t, temp);
                }
            }

            ifft_mfa_truncate_sqrt2(Ccoeffs, n, w, &t1, &t2, &t3, sqrt, trunc);
        }
        else
        {
            for (l = 0; l < trunc; l++)
            {
                _dot(Ccoeffs[l], Acoeffs + (i*K + 0)*4*n + l, 4*n,
                                 Bcoeffs + (0*N + j)*4*n + l, N*4*n,
                     K, limbs, t, temp);
            }

            ifft_truncate_sqrt2(Ccoeffs, n, w, &t1, &t2, &t3, trunc);
        }

        /* pointwise scalar division by 4*n */
        for (l = 0; l < trunc; l++)
        {
            mpn_div_2expmod_2expp1(Ccoeffs[l], Ccoeffs[l], limbs, depth + 2);
            mpn_normmod_2expp1(Ccoeffs[l], limbs);
        }

        fft_combine_bits_fmpz(fmpz_mat_entry(C, i, j), Ccoeffs,
                                      j1 + j2 - 1, bits1, limbs, Climbs, sign);
    }

    flint_free(temp);
    flint_free(coeffs);
}

/*
For nonnegative input:
    just like mpn_mul_main but
        bits = (n*w - (depth + 1))/2;
    is replaced by
        bits = (n*w - (depth + 1 + FLINT_CLOG2(k)))/2;

    the inputs are viewed a polynomials evaluated at 2^bits, with
        all coeffs in [0, 2^bits)

For general input:
    the inputs are viewed a polynomials evaluated at 2^bits, either with
        all coeffs in [0, 2^bits), or
        all coeffs in (-2^bits, 0]

    (Balancing as (-2^(bits-1), 2^(bits-1)] hardly seems worth the effort.)

In either case, the product of the matrices of these polynomials has
coefficients bounded in absolute value by

    k*2^(depth+1)*2^(2*bits)

We want this to be less than 2^(n*w) in the unsigned case and less than
2^(n*w)/2 in the signed case.
*/
void _fmpz_mat_mul_fft(
    fmpz_mat_t C,
    const fmpz_mat_t A, slong abits,
    const fmpz_mat_t B, slong bbits,
    int sign)
{
    slong K = fmpz_mat_ncols(A);
    slong clgK = FLINT_CLOG2(K) + sign;
    slong /*off, */depth = 6;
    slong w = 1;
    slong n = WORD(1) << depth;
    ulong bits = (n*w - (depth + 1 + clgK))/2;
    ulong bits1 = FLINT_MAX(abits, WORD(2000));
    ulong bits2 = FLINT_MAX(bbits, WORD(2000));
    /* when viewed as eval of poly at x=2^bits, the max length of the entries */
    slong j1 = (bits1 + bits - 1)/bits;
    /* ditto for B */
    slong j2 = (bits2 + bits - 1)/bits;
    int use_mfa;

    FLINT_ASSERT(sign == 0 || sign == 1);
    FLINT_ASSERT(abits > 0);
    FLINT_ASSERT(bbits > 0);
    FLINT_ASSERT(j1 + j2 - 1 > 2*n);

    /* find initial n, w */
    while (j1 + j2 - 1 > 4*n)
    {
        if (w == 1)
        {
            w = 2;
        }
        else
        {
            depth++;
            w = 1;
            n *= 2;
        }

        bits = (n*w - (depth + 1 + clgK))/2;
        j1 = (bits1 + bits - 1)/bits;
        j2 = (bits2 + bits - 1)/bits;
    }

    if (depth < 11)
    {
        slong wadj = 1;

        /* adjust n and w */
        if (depth < 6)
            wadj = WORD(1) << (6 - depth);

        if (w > wadj)
        {
            /* see if a smaller w will work */
            do {
                w -= wadj;
                bits = (n*w - (depth + 1 + clgK))/2;
                j1 = (bits1 + bits - 1)/bits;
                j2 = (bits2 + bits - 1)/bits;
            } while (j1 + j2 - 1 <= 4*n && w > wadj);
            w += wadj;
        }

        use_mfa = 0;
    }
    else
    {
        use_mfa = 1;
    }

    bits = (n*w - (depth + 1 + clgK))/2;
    j1 = (bits1 + bits - 1)/bits;
    j2 = (bits2 + bits - 1)/bits;

    _fmpz_mat_mul_truncate_sqrt2(C, A, abits, B, bbits, depth, w, j1, j2,
                                                                use_mfa, sign);
}

void fmpz_mat_mul_fft(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc;
    slong abits, bbits;
    int sign;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    if (ar == 0 || br == 0 || bc == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    abits = fmpz_mat_max_bits(A);
    bbits = fmpz_mat_max_bits(B);

    sign = 0;
    if (abits < 0)
    {
        sign = 1;
        abits = -abits;
    }

    if (bbits < 0)
    {
        sign = 1;
        bbits = -bbits;
    }

    if (abits == 0 || bbits == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    _fmpz_mat_mul_fft(C, A, abits, B, bbits, sign);
}

