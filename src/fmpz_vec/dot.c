/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_extras.h"
#include "fmpz_vec.h"

void
_fmpz_vec_dot_general_naive(fmpz_t res, const fmpz_t initial,
    int subtract, const fmpz * a, const fmpz * b, int reverse, slong len)
{
    slong i;

    if (initial == NULL)
    {
        if (len == 0)
        {
            fmpz_zero(res);
            return;
        }

        fmpz_mul(res, a, reverse ? b + len - 1 : b);

        if (subtract)
        {
            fmpz_neg(res, res);
            for (i = 1; i < len; i++)
                fmpz_submul(res, a + i, reverse ? b + len - 1 - i : b + i);
        }
        else
        {
            for (i = 1; i < len; i++)
                fmpz_addmul(res, a + i, reverse ? b + len - 1 - i : b + i);
        }
    }
    else
    {
        if (res != initial)
            fmpz_set(res, initial);

        if (subtract)
            for(i = 0; i < len; i++)
                fmpz_submul(res, a + i, reverse ? b + len - 1 - i : b + i);
        else
            for(i = 0; i < len; i++)
                fmpz_addmul(res, a + i, reverse ? b + len - 1 - i : b + i);
    }
}

#define INITIAL_ALLOC 64

#define FMPZ_GET_MPN(xptr, xn, xneg, xtmplimb, x) \
    if (!COEFF_IS_MPZ(x)) \
    { \
        (xtmplimb) = FLINT_ABS(x); \
        (xptr) = &(xtmplimb); \
        (xn) = 1; \
        (xneg) = ((x) < 0); \
    } \
    else \
    { \
        mpz_ptr __z1 = COEFF_TO_PTR(x); \
        (xptr) = __z1->_mp_d; \
        (xn) = FLINT_ABS(__z1->_mp_size); \
        (xneg) = (__z1->_mp_size < 0); \
    }

/* (s,sn) = (a,an) + (b,bn). Allows an == 0 but not bn == 0. */
#define MPN_ADD(s, sn, a, an, b, bn) \
    do { \
        if ((an) == 0) \
        { \
            FLINT_SWAP(nn_ptr, s, b); \
            (sn) = (bn); \
        } \
        else if ((an) >= (bn)) \
        { \
            ulong __cy; \
            (s)[(an)] = __cy = mpn_add((s), (a), (an), (b), (bn)); \
            (sn) = (an) + (__cy != 0); \
        } \
        else \
        { \
            ulong __cy; \
            (s)[(bn)] = __cy = mpn_add((s), (b), (bn), (a), (an)); \
            (sn) = (bn) + (__cy != 0); \
        } \
    } while (0)

/* (s,sn) = (s,sn) + (a,an) * b. Allows sn == 0 but not an == 0. */
#define MPN_ADDMUL_1(s, sn, a, an, b) \
    do { \
        ulong __cy; \
        if ((sn) >= (an)) \
        { \
            FLINT_ASSERT((an) != 0); \
            __cy = mpn_addmul_1((s), (a), (an), (b)); \
            if ((sn) > (an)) \
                __cy = mpn_add_1((s) + (an), (s) + (an), (sn) - (an), __cy); \
            (s)[(sn)] = __cy; \
            (sn) += (__cy != 0); \
        } \
        else \
        { \
            (s)[(an)] = mpn_mul_1((s) + (sn), (a) + (sn), (an) - (sn), (b)); \
            if ((sn) != 0) \
            { \
                __cy = mpn_addmul_1((s), (a), (sn), (b)); \
                (s)[(an)] += mpn_add_1((s) + (sn), (s) + (sn), (an) - (sn), __cy); \
            } \
            (sn) = (an) + ((s)[(an)] != 0); \
        } \
    } while (0) \


FLINT_STATIC_NOINLINE
void _fmpz_set_mpn(fmpz_t res, nn_srcptr x, slong xn, int neg)
{
    if (xn <= 1 && x[0] <= COEFF_MAX)
    {
        _fmpz_demote(res);
        *res = neg ? -x[0] : x[0];
    }
    else
    {
        fmpz_set_mpn_large(res, x, xn, neg);
    }
}

void
_fmpz_vec_dot_general(fmpz_t res, const fmpz_t initial, int subtract,
            const fmpz * a, const fmpz * b, int reverse, slong len)
{
    ulong tmp1[INITIAL_ALLOC + 2];
    ulong tmp2[INITIAL_ALLOC + 2];
    ulong tmp3[INITIAL_ALLOC + 2];
    slong alloc = INITIAL_ALLOC;
    slong new_alloc;

    /* We maintain separate sums for small terms, large positive terms,
       and large negative terms, the idea being to have fewer
       adjustments in the main loop in exchange for some added
       complexity combining things in the end. Should profile
       alternative strategies. */
    ulong s0 = 0, s1 = 0, s2 = 0;
    nn_ptr neg = tmp1;
    nn_ptr pos = tmp2;
    slong posn = 0, negn = 0;

    /* Temporary space for products. */
    nn_ptr t = tmp3;
    slong tn;

    nn_ptr tmp_heap = NULL;

    slong i;

    if (len <= 1)
    {
        if (initial == NULL)
        {
            if (len == 1)
            {
                fmpz_mul(res, a, b);
                if (subtract)
                    fmpz_neg(res, res);
            }
            else
                fmpz_zero(res);
        }
        else
        {
            if (res != initial)
                fmpz_set(res, initial);

            if (subtract)
            {
                if (len == 1)
                    fmpz_submul(res, a, b);
            }
            else
            {
                if (len == 1)
                    fmpz_addmul(res, a, b);
            }
        }
        return;
    }

    if (initial != NULL)
    {
        fmpz ca;
        ulong atmp;
        nn_srcptr ap;
        slong an;
        int aneg;

        ca = *initial;
        FMPZ_GET_MPN(ap, an, aneg, atmp, ca);

        if (an <= 2)
        {
            s0 = ap[0];
            if (an == 2)
                s1 = ap[1];

            if (aneg ^ subtract)
                sub_dddmmmsss(s2, s1, s0, 0, 0, 0, 0, s1, s0);
        }
        else
        {
            if (an > INITIAL_ALLOC)
            {
                new_alloc = an + 4;

                tmp_heap = flint_malloc(3 * (new_alloc + 2) * sizeof(ulong));

                t = tmp_heap;
                pos = t + (new_alloc + 2);
                neg = pos + (new_alloc + 2);

                alloc = new_alloc;
            }

            if (aneg ^ subtract)
            {
                flint_mpn_copyi(neg, ap, an);
                negn = an;
            }
            else
            {
                flint_mpn_copyi(pos, ap, an);
                posn = an;
            }
        }
    }

    for (i = 0; i < len; i++)
    {
        fmpz ca, cb;
        ulong atmp, btmp;
        nn_srcptr ap, bp;
        slong an, bn;
        ulong cy;
        int aneg, bneg;

        ca = a[i];
        if (ca == 0)
            continue;

        cb = reverse ? b[len - 1 - i] : b[i];
        if (cb == 0)
            continue;

        if (!COEFF_IS_MPZ(ca) && !COEFF_IS_MPZ(cb))
        {
            ulong hi, lo;
            smul_ppmm(hi, lo, ca, cb);
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, FLINT_SIGN_EXT(hi), hi, lo);
            continue;
        }

        FMPZ_GET_MPN(ap, an, aneg, atmp, ca);
        FMPZ_GET_MPN(bp, bn, bneg, btmp, cb);
        tn = an + bn;

        if (tn > alloc)
        {
            nn_ptr p1, p2, p3;

            new_alloc = FLINT_MAX(3 * alloc / 2, tn + 4);

            p1 = flint_malloc(3 * (new_alloc + 2) * sizeof(ulong));
            p2 = p1 + (new_alloc + 2);
            p3 = p2 + (new_alloc + 2);

            flint_mpn_copyi(p2, pos, posn);
            flint_mpn_copyi(p3, neg, negn);
            t = p1;
            pos = p2;
            neg = p3;

            FLINT_SWAP(nn_ptr, tmp_heap, p1);

            if (p1 != NULL)
                flint_free(p1);

            alloc = new_alloc;
        }

        if (an < bn)
        {
            FLINT_SWAP(nn_srcptr, ap, bp);
            FLINT_SWAP(slong, an, bn);
        }

        if (bn == 1)
        {
            ulong b0 = bp[0];

            if (aneg ^ bneg)
                MPN_ADDMUL_1(neg, negn, ap, an, b0);
            else
                MPN_ADDMUL_1(pos, posn, ap, an, b0);
            continue;
        }

        if (ap == bp && an == bn)
        {
            flint_mpn_sqr(t, ap, an);
            cy = t[tn - 1];
        }
        else
        {
            cy = flint_mpn_mul(t, ap, an, bp, bn);
        }

        tn -= (cy == 0);

        if (aneg ^ bneg)
            MPN_ADD(neg, negn, neg, negn, t, tn);
        else
            MPN_ADD(pos, posn, pos, posn, t, tn);
    }

    /* There are only small terms. */
    if (posn == 0 && negn == 0)
    {
        if (subtract)
            sub_dddmmmsss(s2, s1, s0, 0, 0, 0, s2, s1, s0);

        fmpz_set_signed_uiuiui(res, s2, s1, s0);
        return;
    }

    /* Add small terms to large ones. */
    if (((slong) s2 >= WORD(0)))
    {
        t[2] = s2;
        t[1] = s1;
        t[0] = s0;
        MPN_ADD(pos, posn, pos, posn, t, 3);
    }
    else
    {
        sub_dddmmmsss(t[2], t[1], t[0], 0, 0, 0, s2, s1, s0);
        MPN_ADD(neg, negn, neg, negn, t, 3);
    }

    MPN_NORM(pos, posn);
    MPN_NORM(neg, negn);

    if (negn == 0)
    {
        _fmpz_set_mpn(res, pos, posn, 0 ^ subtract);
    }
    else if (posn == 0)
    {
        _fmpz_set_mpn(res, neg, negn, 1 ^ subtract);
    }
    else
    {
        /* Do subtraction */
        int tneg = 0;

        if (posn > negn)
        {
            tn = posn;
        }
        else if (negn > posn)
        {
            tn = negn;
            tneg = 1;
        }
        else if (posn != 0)
        {
            tn = posn;
            if (mpn_cmp(pos, neg, tn) < 0)
                tneg = 1;
        }

        if (tneg)
            mpn_sub(t, neg, negn, pos, posn);
        else
            mpn_sub(t, pos, posn, neg, negn);

        MPN_NORM(t, tn);
        _fmpz_set_mpn(res, t, tn, tneg ^ subtract);
    }

    if (tmp_heap != NULL)
        flint_free(tmp_heap);
}
