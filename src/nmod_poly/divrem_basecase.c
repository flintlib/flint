/*
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly/impl.h"

/* addmul_fits: a + b*c fits in 16/32/64 bits */
/* addmul_addmul_fits: a + b*c + d*e fits in 16/32/64 bits */

#if FLINT_BITS == 64

#define NMOD_ADDMUL_FITS_HALFWORD(mod) ((mod.n) <= UWORD(65535))
#define NMOD_ADDMUL_ADDMUL_FITS_HALFWORD(mod) ((mod.n) <= UWORD(46341))
#define NMOD_ADDMUL_FITS_WORD(mod) ((mod.n) <= UWORD(4294967295))
#define NMOD_ADDMUL_ADDMUL_FITS_WORD(mod) ((mod.n) <= UWORD(3037000500))

#else

#define NMOD_ADDMUL_FITS_HALFWORD(mod) ((mod.n) <= UWORD(255))
#define NMOD_ADDMUL_ADDMUL_FITS_HALFWORD(mod) ((mod.n) <= UWORD(181))
#define NMOD_ADDMUL_FITS_WORD(mod) ((mod.n) <= UWORD(65535))
#define NMOD_ADDMUL_ADDMUL_FITS_WORD(mod) ((mod.n) <= UWORD(46341))

#endif

/* Lemire, Kaser & Kurz */
static ulong _mod1_preinvert(ulong n)
{
    return UWORD_MAX / n + 1;
}

/* Shoup */
static ulong _mod2_preinvert(ulong n)
{
    return n_mulmod_precomp_shoup(1, n);
}

static ulong _mod1(ulong x, ulong n, ulong ninv)
{
    ulong l = ninv * x;

#if FLINT_BITS != 64 || !defined(__GNUC__)
    ulong hi, lo;
    umul_ppmm(hi, lo, l, n);
    return hi;
#else
    return ((__uint128_t) l * n ) >> 64;
#endif
}

static ulong _mulmod1(ulong a, ulong b, ulong n, ulong ninv)
{
    return _mod1(a * b, n, ninv);
}

static ulong _mod2_fast(ulong x, ulong n, ulong ninv)
{
    return x - n_mulhi(ninv, x) * n;
}

static ulong _mod2(ulong x, ulong n, ulong ninv)
{
    ulong y = _mod2_fast(x, n, ninv);
    if (y >= n)
        y -= n;
    return y;
}

static ulong _mulmod2(ulong a, ulong b, ulong n, ulong ninv)
{
    return _mod2(a * b, n, ninv);
}

void _nmod_poly_divrem_q0_preinv1(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, nn_srcptr B, slong lenA, ulong invL, nmod_t mod)
{
    if (lenA == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        Q[0] = nmod_mul(A[lenA-1], invL, mod);

        if (R == A)
        {
            _nmod_vec_scalar_addmul_nmod(R, B, lenA - 1, nmod_neg(Q[0], mod), mod);
        }
        else
        {
            _nmod_vec_scalar_mul_nmod(R, B, lenA - 1, Q[0], mod);
            _nmod_vec_sub(R, A, R, lenA - 1, mod);
        }
    }
}

/* For GCC and Zen 3, the performance of the following critical function is
   extremely finicky, with an immediate 10% slowdown if one reorders some
   instruction in the main loop, makes it inline, etc. Just try to profile
   and do whatever works.
*/

FLINT_STATIC_NOINLINE
void _nmod_poly_divrem_q1_preinv1_fullword(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          ulong invL, nmod_t mod)
{
    slong i;
    ulong t, q0, q1, t1, t0, s1, s0;

    FLINT_ASSERT(lenA == lenB + 1);

    q1 = _nmod_mul_fullword(A[lenA-1], invL, mod);
    t  = _nmod_mul_fullword(q1, B[lenB-2], mod);
    t  = nmod_sub(t, A[lenA-2], mod);
    q0 = _nmod_mul_fullword(t, invL, mod);
    Q[0] = nmod_neg(q0, mod);
    Q[1] = q1;
    q1 = nmod_neg(q1, mod);
    /* R = A + (q1*x + q0)*B */

    /* Hack: nmod_addmul is faster than _nmod_mull_fullword + nmod_add */
    R[0] = nmod_addmul(A[0], q0, B[0], mod);

    /* We want to compute r = a + q0*b + q1*c where a, b, c, q0, q1 <= n - 1.
       If (1 + q0 + q1) < 2^FLINT_BITS, then r < 2^FLINT_BITS * (n-1).
       In this case, r fits in two limbs without overflow and the high limb is
       already reduced mod n. This will almost always be the case when n is
       close to 2^(FLINT_BITS-1), and happens with decent probability when
       n is close to 2^FLINT_BITS too. */
    if (q0 + q1 + 1 >= q0)  /* Checks that q0 + (q1 + 1) doesn't wrap around. */
    {
        for (i = 1; i < lenB - 1; i++)
        {
            umul_ppmm(t1, t0, q1, B[i - 1]);
            add_ssaaaa(t1, t0, t1, t0, 0, A[i]);
            umul_ppmm(s1, s0, q0, B[i]);
            add_ssaaaa(t1, t0, t1, t0, s1, s0);
            FLINT_ASSERT(t1 < mod.n);
            NMOD_RED2_FULLWORD(R[i], t1, t0, mod);
        }
    }
    else
    {
        for (i = 1; i < lenB - 1; i++)
        {
            umul_ppmm(t1, t0, q1, B[i - 1]);
            add_ssaaaa(t1, t0, t1, t0, 0, A[i]);
            umul_ppmm(s1, s0, q0, B[i]);
            add_ssaaaa(t1, t0, t1, t0, 0, s0);
            t1 = nmod_add(t1, s1, mod);
            FLINT_ASSERT(t1 < mod.n);
            NMOD_RED2_FULLWORD(R[i], t1, t0, mod);
        }
    }
}

void _nmod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          ulong invL, nmod_t mod)
{
    slong i;
    ulong t, q0, q1, t1, t0, s1, s0;

    FLINT_ASSERT(lenA == lenB + 1);

    if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
        return;
    }

    if (NMOD_ADDMUL_FITS_HALFWORD(mod))
    {
        ulong n = mod.n;
        ulong ninv = _mod1_preinvert(n);

        q1 = _mulmod1(A[lenA-1], invL, n, ninv);
        t  = _mulmod1(q1, B[lenB-2], n, ninv);
        t  = nmod_sub(t, A[lenA-2], mod);
        q0 = _mulmod1(t, invL, n, ninv);
        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);
        /* R = A + (q1*x + q0)*B */

        R[0] = _mod1(A[0] + q0 * B[0], n, ninv);

        if (NMOD_ADDMUL_ADDMUL_FITS_HALFWORD(mod))
        {
            for (i = 1; i < lenB - 1; i++)
                R[i] = _mod1(A[i] + q1*B[i - 1] + q0*B[i], n, ninv);
        }
        else
        {
            for (i = 1; i < lenB - 1; i++)
                R[i] = _mod1(_mod1(A[i] + q1*B[i - 1], n, ninv) + q0*B[i], n, ninv);
        }
    }
    else if (NMOD_ADDMUL_FITS_WORD(mod))
    {
        ulong n = mod.n;
        ulong ninv = _mod2_preinvert(n);

        q1 = _mulmod2(A[lenA-1], invL, n, ninv);
        t  = _mulmod2(q1, B[lenB-2], n, ninv);
        t  = nmod_sub(t, A[lenA-2], mod);
        q0 = _mulmod2(t, invL, n, ninv);
        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);
        /* R = A + (q1*x + q0)*B */

        R[0] = _mod2(A[0] + q0 * B[0], n, ninv);

        /* In the second branch, the _mod2_fast maps
               [0, ..., (n-1) + (n-1)^2] -> [0, 2n-1]
           and the _mod2 maps
                [0, 2n-1 + (n-1)^2] -> [0,n-1]. */

        if (NMOD_ADDMUL_ADDMUL_FITS_WORD(mod))
            for (i = 1; i < lenB - 1; i++)
                R[i] = _mod2(A[i] + q1*B[i - 1] + q0*B[i], n, ninv);
        else
            for (i = 1; i < lenB - 1; i++)
                R[i] = _mod2(_mod2_fast(A[i] + q1*B[i - 1], n, ninv) + q0*B[i], n, ninv);
    }
    else if (NMOD_BITS(mod) != FLINT_BITS)
    {
        FLINT_ASSERT(lenA == lenB + 1);

        q1 = nmod_mul(A[lenA-1], invL, mod);
        t  = nmod_mul(q1, B[lenB-2], mod);
        t  = nmod_sub(t, A[lenA-2], mod);
        q0 = nmod_mul(t, invL, mod);
        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);
        /* R = A + (q1*x + q0)*B */

        R[0] = nmod_addmul(A[0], q0, B[0], mod);

        if (NMOD_BITS(mod) != FLINT_BITS)
        {
            for (i = 1; i < lenB - 1; i++)
            {
                umul_ppmm(t1, t0, q1, B[i - 1]);
                umul_ppmm(s1, s0, q0, B[i]);
                add_ssaaaa(t1, t0, t1, t0, 0, A[i]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                t1 = FLINT_MIN(t1, t1 - mod.n);
                FLINT_ASSERT(t1 < mod.n);
                NMOD_RED2(R[i], t1, t0, mod);
            }
        }
    }
    else
    {
        _nmod_poly_divrem_q1_preinv1_fullword(Q, R, A, lenA, B, lenB, invL, mod);
    }
}

static void
_nmod_poly_divrem_basecase_preinv1_1(nn_ptr Q, nn_ptr R, nn_ptr W,
                             nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                             ulong invL,
                             nmod_t mod)
{
    slong iR;
    nn_ptr ptrQ = Q - lenB + 1;
    nn_ptr R1 = W;

    flint_mpn_copyi(R1, A, lenA);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (R1[iR] == 0)
        {
            ptrQ[iR] = WORD(0);
        }
        else
        {
            ptrQ[iR] = n_mulmod2_preinv(R1[iR], invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const ulong c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R1 + iR - lenB + 1, B, lenB - 1, c);
            }
        }
    }

    if (lenB > 1)
        _nmod_vec_reduce(R, R1, lenB - 1, mod);
}

static void
_nmod_poly_divrem_basecase_preinv1_2(nn_ptr Q, nn_ptr R, nn_ptr W,
                             nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                             ulong invL,
                             nmod_t mod)
{
    slong iR, i;
    nn_ptr B2 = W, R2 = W + 2*(lenB - 1), ptrQ = Q - lenB + 1;

    for (i = 0; i < lenB - 1; i++)
    {
        B2[2 * i] = B[i];
        B2[2 * i + 1] = 0;
    }
    for (i = 0; i < lenA; i++)
    {
        R2[2 * i] = A[i];
        R2[2 * i + 1] = 0;
    }

    for (iR = lenA - 1; iR >= lenB - 1; )
    {
        ulong r =
            n_ll_mod_preinv(R2[2 * iR + 1], R2[2 * iR], mod.n, mod.ninv);

        while ((iR + 1 >= lenB) && (r == WORD(0)))
        {
            ptrQ[iR--] = WORD(0);
            if (iR + 1 >= lenB)
                r = n_ll_mod_preinv(R2[2 * iR + 1], R2[2 * iR], mod.n,
                                    mod.ninv);
        }

        if (iR + 1 >= lenB)
        {
            ptrQ[iR] = n_mulmod2_preinv(r, invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const ulong c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R2 + 2 * (iR - lenB + 1), B2, 2 * lenB - 2, c);
            }
            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        R[iR] = n_ll_mod_preinv(R2[2*iR+1], R2[2*iR], mod.n, mod.ninv);
}

static void
_nmod_poly_divrem_basecase_preinv1_3(nn_ptr Q, nn_ptr R, nn_ptr W,
                                     nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                                     ulong invL,
                                     nmod_t mod)
{
    slong iR, i;
    nn_ptr B3 = W, R3 = W + 3*(lenB - 1), ptrQ = Q - lenB + 1;

    for (i = 0; i < lenB - 1; i++)
    {
        B3[3 * i] = B[i];
        B3[3 * i + 1] = 0;
        B3[3 * i + 2] = 0;
    }
    for (i = 0; i < lenA; i++)
    {
        R3[3 * i] = A[i];
        R3[3 * i + 1] = 0;
        R3[3 * i + 2] = 0;
    }

    for (iR = lenA - 1; iR >= lenB - 1; )
    {
        ulong r =
            n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                             R3[3 * iR], mod.n, mod.ninv);

        while ((iR + 1 >= lenB) && (r == WORD(0)))
        {
            ptrQ[iR--] = WORD(0);
            if (iR + 1 >= lenB)
                r = n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                                     R3[3 * iR], mod.n, mod.ninv);
        }

        if (iR + 1 >= lenB)
        {
            ptrQ[iR] = n_mulmod2_preinv(r, invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const ulong c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R3 + 3 * (iR - lenB + 1), B3, 3 * lenB - 3, c);
            }
            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        R[iR] = n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                                 R3[3 * iR], mod.n, mod.ninv);
}

static
slong NMOD_DIVREM_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const flint_bitcnt_t bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);

    if (bits <= FLINT_BITS)
        return lenA;
    else if (bits <= 2 * FLINT_BITS)
        return 2*(lenA + lenB - 1);
    else
        return 3*(lenA + lenB - 1);
}

void
_nmod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R,
                           nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                           ulong invB,
                           nmod_t mod)
{
    if (lenA == lenB + 1)
    {
        _nmod_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invB, mod);
    }
    else if (lenA == lenB)
    {
        _nmod_poly_divrem_q0_preinv1(Q, R, A, B, lenB, invB, mod);
    }
    else if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invB, mod);
    }
    else
    {
        nn_ptr W;
        TMP_INIT;
        slong bits = 2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);

        TMP_START;
        W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod)*sizeof(ulong));

        if (bits <= FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_1(Q, R, W, A, lenA, B, lenB, invB, mod);
        else if (bits <= 2 * FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_2(Q, R, W, A, lenA, B, lenB, invB, mod);
        else
            _nmod_poly_divrem_basecase_preinv1_3(Q, R, W, A, lenA, B, lenB, invB, mod);

        TMP_END;
    }
}

void
_nmod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA,
                                  nn_srcptr B, slong lenB, nmod_t mod)
{
    ulong invB;

    invB = (B[lenB - 1] == 1) ? 1 : n_invmod(B[lenB - 1], mod.n);
    _nmod_poly_divrem_basecase_preinv1(Q, R, A, lenA, B, lenB, invB, mod);
}

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tQ, tR;
    nn_ptr q, r;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            nmod_poly_zero(R);
            return;
        } else
        {
            flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_divrem). Division by zero.");
        }
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(tQ, A->mod.n, A->mod.ninv, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == B)
    {
        nmod_poly_init2_preinv(tR, B->mod.n, B->mod.ninv, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, tQ);
        nmod_poly_clear(tQ);
    }
    if (R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }

    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
