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
        ulong ninv = n_lemire_precomp_preinv(n, mod.ninv);

        q1 = n_mod_lemire(A[lenA-1] * invL, n, ninv);
        t  = n_mod_lemire(q1 * B[lenB-2], n, ninv);
        t  = nmod_sub(t, A[lenA-2], mod);
        q0 = n_mod_lemire(t * invL, n, ninv);
        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);
        /* R = A + (q1*x + q0)*B */

        R[0] = n_mod_lemire(A[0] + q0 * B[0], n, ninv);

        if (NMOD_ADDMUL_ADDMUL_FITS_HALFWORD(mod))
        {
            for (i = 1; i < lenB - 1; i++)
                R[i] = n_mod_lemire(A[i] + q1*B[i - 1] + q0*B[i], n, ninv);
        }
        else
        {
            for (i = 1; i < lenB - 1; i++)
                R[i] = n_mod_lemire(n_mod_lemire(A[i] + q1*B[i - 1], n, ninv) + q0*B[i], n, ninv);
        }
    }
    else if (NMOD_ADDMUL_FITS_WORD(mod))
    {
        ulong n = mod.n;
        ulong ninv = n_barrett_precomp_preinv(n, mod.ninv);

        q1 = n_mod_barrett(A[lenA-1] * invL, n, ninv);
        t  = n_mod_barrett(q1 * B[lenB-2], n, ninv);
        t  = nmod_sub(t, A[lenA-2], mod);
        q0 = n_mod_barrett(t * invL, n, ninv);
        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);
        /* R = A + (q1*x + q0)*B */

        R[0] = n_mod_barrett(A[0] + q0 * B[0], n, ninv);

        /* In the second branch, the lazy mod maps
               [0, ..., (n-1) + (n-1)^2] -> [0, 2n-1]
           and the _mod2 maps
                [0, 2n-1 + (n-1)^2] -> [0,n-1]. */

        if (NMOD_ADDMUL_ADDMUL_FITS_WORD(mod))
            for (i = 1; i < lenB - 1; i++)
                R[i] = n_mod_barrett(A[i] + q1*B[i - 1] + q0*B[i], n, ninv);
        else
            for (i = 1; i < lenB - 1; i++)
                R[i] = n_mod_barrett(n_mod_barrett_lazy(A[i] + q1*B[i - 1], n, ninv) + q0*B[i], n, ninv);
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

        if (NMOD_CAN_USE_FMMA_BOUNDED(mod))
        {
            /* a residue plus two products, below 2 n^2 < 2^(2 nu + 1)
               <= 2^(nu + FLINT_BITS - 2) as for _nmod_fmma_bounded */
            for (i = 1; i < lenB - 1; i++)
            {
                umul_ppmm(t1, t0, q1, B[i - 1]);
                umul_ppmm(s1, s0, q0, B[i]);
                add_ssaaaa(t1, t0, t1, t0, 0, A[i]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                NMOD_RED2_BOUNDED(R[i], t1, t0, mod);
            }
        }
        else
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
_nmod_poly_divrem_basecase_preinv1_1(nn_ptr Q, nn_ptr R,
                             nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                             ulong invL,
                             nmod_t mod)
{
    slong iR, i;
    nn_ptr R1, ptrQ = Q - lenB + 1;
    ulong r, c;
    int unreduced_fits_halflimb;
    int reduced_fits_quarterlimb;

    FLINT_ASSERT(NMOD_BITS(mod) <= FLINT_BITS / 2);

    ulong unreduced_bound = (mod.n - 1) * (mod.n - 1) * (lenA - lenB + 1);

    unreduced_fits_halflimb = unreduced_bound < UWORD(1) << (FLINT_BITS / 2);
    reduced_fits_quarterlimb = NMOD_BITS(mod) <= FLINT_BITS / 4;

    ulong npre = n_barrett_precomp_preinv(mod.n, mod.ninv);
    ulong npre2 = n_lemire_precomp_preinv(mod.n, mod.ninv);

    TMP_INIT;
    TMP_START;
    R1 = TMP_ALLOC(lenA * sizeof(ulong));

    flint_mpn_copyi(R1, A, lenA);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (unreduced_fits_halflimb)
            r = n_mod_lemire(R1[iR], mod.n, npre2);
        else
            r = n_mod_barrett(R1[iR], mod.n, npre);

        if (r == 0)
        {
            ptrQ[iR] = 0;
        }
        else
        {
            if (invL == 1)
                ptrQ[iR] = r;
            else if (reduced_fits_quarterlimb)
                ptrQ[iR] = n_mod_lemire(invL * r, mod.n, npre2);
            else
                ptrQ[iR] = n_mod_barrett(invL * r, mod.n, npre);

            if (lenB > 1)
            {
                c = mod.n - ptrQ[iR];
                _nmod_vec_nored_scalar_addmul_halflimb(R1 + iR - lenB + 1, B, lenB - 1, c);
            }
        }
    }

    if (lenB > 1)
    {
        if (unreduced_fits_halflimb)
        {
            for (i = 0; i < lenB - 1; i++)
                R[i] = n_mod_lemire(R1[i], mod.n, npre2);
        }
        else
        {
            for (i = 0; i < lenB - 1; i++)
                R[i] = n_mod_barrett(R1[i], mod.n, npre);
        }
    }

    TMP_END;
}

/* Note: we can do better than n_ll_mod_preinv when the high limb
   of the double-limb sum is not reduced, but this benefits only a narrow
   set of parameters, so it is not currently implemented. */
static void
_nmod_poly_divrem_basecase_preinv1_2(nn_ptr Q, nn_ptr R,
                             nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                             ulong invL,
                             nmod_t mod)
{
    slong iR, i;
    nn_ptr R2, ptrQ = Q - lenB + 1;
    ulong r, c;
    int halflimb, bounded;

    TMP_INIT;
    TMP_START;
    R2 = TMP_ALLOC(2 * lenA * sizeof(ulong));

    for (i = 0; i < lenA; i++)
    {
        R2[2 * i] = A[i];
        R2[2 * i + 1] = 0;
    }

    halflimb = NMOD_IS_HALFLIMB(mod);

    /* each entry is a residue plus at most lenA - lenB + 1 products, below
       2^(2 nu + bitcount(lenA - lenB + 1)): at most 2^(nu + FLINT_BITS - 2)
       if nu + bitcount(lenA - lenB + 1) <= FLINT_BITS - 2, with nu >= 3 */
    bounded = (NMOD_BITS(mod) >= 3) && (NMOD_BITS(mod)
               + FLINT_BIT_COUNT(lenA - lenB + 1) <= FLINT_BITS - 2);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (bounded)
            NMOD_RED2_BOUNDED(r, R2[2 * iR + 1], R2[2 * iR], mod);
        else
            r = n_ll_mod_preinv(R2[2 * iR + 1], R2[2 * iR], mod.n, mod.ninv);

        if (r == 0)
        {
            ptrQ[iR] = 0;
        }
        else
        {
            if (invL == 1)
                ptrQ[iR] = r;
            else
                ptrQ[iR] = nmod_mul(r, invL, mod);

            if (lenB > 1)
            {
                c = mod.n - ptrQ[iR];

                if (halflimb)
                    _nmod_vec_nored_ll_scalar_addmul_halflimb(R2 + 2 * (iR - lenB + 1), B, lenB - 1, c);
                else
                    _nmod_vec_nored_ll_scalar_addmul(R2 + 2 * (iR - lenB + 1), B, lenB - 1, c);
            }
        }
    }

    if (bounded)
        for (iR = 0; iR < lenB - 1; iR++)
            NMOD_RED2_BOUNDED(R[iR], R2[2*iR+1], R2[2*iR], mod);
    else
        for (iR = 0; iR < lenB - 1; iR++)
            R[iR] = n_ll_mod_preinv(R2[2*iR+1], R2[2*iR], mod.n, mod.ninv);

    TMP_END;
}

static void
_nmod_poly_divrem_basecase_preinv1_3(nn_ptr Q, nn_ptr R,
                                     nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                                     ulong invL,
                                     nmod_t mod)
{
    slong iR, i;
    nn_ptr R3, ptrQ = Q - lenB + 1;
    ulong r, c;

    TMP_INIT;
    TMP_START;
    R3 = TMP_ALLOC(3 * lenA * sizeof(ulong));

    for (i = 0; i < lenA; i++)
    {
        R3[3 * i] = A[i];
        R3[3 * i + 1] = 0;
        R3[3 * i + 2] = 0;
    }

    /* each entry is a coefficient of A plus at most lenA - lenB + 1
       products of residues, below (lenA - lenB + 1) B^2 */
    nmod_red3_precomp_struct R3P;
    _nmod_red3_precomp_init(&R3P, lenA - lenB + 1, mod);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (R3P.nosub)
            NMOD_RED3_FOLD_FULLWORD_NOSUB(r, R3[3 * iR + 2], R3[3 * iR + 1],
                                          R3[3 * iR], R3P.b1, R3P.b2, mod);
        else
            r = _nmod_red3_precomp(R3[3 * iR + 2], R3[3 * iR + 1], R3[3 * iR], &R3P, mod);

        if (r == 0)
        {
            ptrQ[iR] = 0;
        }
        else
        {
            if (invL == 1)
                ptrQ[iR] = r;
            else
                ptrQ[iR] = nmod_mul(r, invL, mod);

            if (lenB > 1)
            {
                c = mod.n - ptrQ[iR];
                _nmod_vec_nored_lll_scalar_addmul(R3 + 3 * (iR - lenB + 1), B, lenB - 1, c);
            }
        }
    }

    if (R3P.nosub)
        for (iR = 0; iR < lenB - 1; iR++)
            NMOD_RED3_FOLD_FULLWORD_NOSUB(R[iR], R3[3 * iR + 2], R3[3 * iR + 1],
                                          R3[3 * iR], R3P.b1, R3P.b2, mod);
    else
        for (iR = 0; iR < lenB - 1; iR++)
            R[iR] = _nmod_red3_precomp(R3[3 * iR + 2], R3[3 * iR + 1], R3[3 * iR], &R3P, mod);

    TMP_END;
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

        slong bits = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(lenA - lenB + 1);

        if (bits <= FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_1(Q, R, A, lenA, B, lenB, invB, mod);
        else if (bits <= 2 * FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_2(Q, R, A, lenA, B, lenB, invB, mod);
        else
            _nmod_poly_divrem_basecase_preinv1_3(Q, R, A, lenA, B, lenB, invB, mod);
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
            flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_divrem_basecase). Division by zero.");
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

    _nmod_poly_divrem_basecase(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

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
