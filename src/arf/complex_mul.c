/*
    Copyright (C) 2014, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"
#include "arf.h"
#include "mpn_extras.h"

int arf_complex_mul_fallback(arf_t e, arf_t f,
        const arf_t a, const arf_t b,
        const arf_t c, const arf_t d,
        slong prec, arf_rnd_t rnd)
{
    int inex1, inex2;

    /* here the operations are ordered to allow aliasing */
    if (arf_is_zero(d))
    {
        /* (a + bi) * c */
        inex2 = arf_mul(f, b, c, prec, rnd);
        inex1 = arf_mul(e, a, c, prec, rnd);
    }
    else if (arf_is_zero(b))
    {
        /* a * (c + di) */
        inex2 = arf_mul(f, a, d, prec, rnd);
        inex1 = arf_mul(e, a, c, prec, rnd);
    }
    else if (arf_is_zero(c))
    {
        /* (a + bi) * di = -bd + adi */
        inex2 = arf_mul(e, a, d, prec, rnd);
        inex1 = arf_neg_mul(f, b, d, prec, rnd);
        arf_swap(e, f);
    }
    else if (arf_is_zero(a))
    {
        /* bi * (c + di) = -bd + bci */
        inex2 = arf_mul(e, b, c, prec, rnd);
        inex1 = arf_neg_mul(f, b, d, prec, rnd);
        arf_swap(e, f);
    }
    else
    {
        arf_t t, u, v, w;

        arf_init(t);
        arf_init(u);
        arf_init(v);
        arf_init(w);

        arf_mul(t, a, c, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul(u, b, d, ARF_PREC_EXACT, ARF_RND_DOWN);

        arf_mul(v, b, c, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul(w, a, d, ARF_PREC_EXACT, ARF_RND_DOWN);

        inex1 = arf_sub(e, t, u, prec, rnd);
        inex2 = arf_add(f, v, w, prec, rnd);

        arf_clear(t);
        arf_clear(u);
        arf_clear(v);
        arf_clear(w);
    }

    return inex1 | (inex2 << 1);
}


int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b,
                                      const arf_t c, const arf_t d,
                                      slong prec, arf_rnd_t rnd)
{
    mp_srcptr ap, bp, cp, dp;
    int asgn, bsgn, csgn, dsgn, inex1, inex2;
    mp_ptr tmp, acp, bdp, adp, bcp;
    mp_size_t an, bn, cn, dn, acn, bdn, adn, bcn, alloc;
    slong shift;
    slong aexp, bexp, cexp, dexp;
    fmpz texp, uexp;

    if (!ARF_IS_LAGOM(a) || !ARF_IS_LAGOM(b) ||
        !ARF_IS_LAGOM(c) || !ARF_IS_LAGOM(d) ||
        ARF_IS_SPECIAL(d) || ARF_IS_SPECIAL(a) ||
        ARF_IS_SPECIAL(b) || ARF_IS_SPECIAL(c))
    {
        return arf_complex_mul_fallback(e, f, a, b, c, d, prec, rnd);
    }

    ARF_GET_MPN_READONLY(ap, an, a);
    asgn = ARF_SGNBIT(a);
    aexp = ARF_EXP(a);

    ARF_GET_MPN_READONLY(bp, bn, b);
    bsgn = ARF_SGNBIT(b);
    bexp = ARF_EXP(b);

    ARF_GET_MPN_READONLY(cp, cn, c);
    csgn = ARF_SGNBIT(c);
    cexp = ARF_EXP(c);

    ARF_GET_MPN_READONLY(dp, dn, d);
    dsgn = ARF_SGNBIT(d);
    dexp = ARF_EXP(d);

    /* Karatsuba multiplication
        e = ac - bd
        f = (a+b)(c+d) - ac - bd */
    if (an >= 20 &&
        cn >= 20 &&
        FLINT_ABS(an - bn) <= 2 &&
        FLINT_ABS(cn - dn) <= 2 &&
        FLINT_ABS(aexp - bexp) <= 64 &&
        FLINT_ABS(cexp - dexp) <= 64)
    {
        fmpz_t za, zb, zc, zd, t, u, v;
        slong abot, bbot, cbot, dbot;

        abot = aexp - an * FLINT_BITS;
        bbot = bexp - bn * FLINT_BITS;
        cbot = cexp - cn * FLINT_BITS;
        dbot = dexp - dn * FLINT_BITS;

        texp = FLINT_MIN(abot, bbot);
        uexp = FLINT_MIN(cbot, dbot);

        fmpz_init(za);
        fmpz_init(zb);
        fmpz_init(zc);
        fmpz_init(zd);
        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(v);

        /* todo: could avoid two copies */
        fmpz_lshift_mpn(za, ap, an, asgn, abot - texp);
        fmpz_lshift_mpn(zb, bp, bn, bsgn, bbot - texp);
        fmpz_lshift_mpn(zc, cp, cn, csgn, cbot - uexp);
        fmpz_lshift_mpn(zd, dp, dn, dsgn, dbot - uexp);

        fmpz_add(t, za, zb);
        fmpz_add(v, zc, zd);
        fmpz_mul(u, t, v);
        fmpz_mul(t, za, zc);
        fmpz_mul(v, zb, zd);
        fmpz_sub(u, u, t);
        fmpz_sub(u, u, v);
        fmpz_sub(t, t, v);

        texp += uexp;
        inex1 = arf_set_round_fmpz_2exp(e, t, &texp, prec, rnd);
        inex2 = arf_set_round_fmpz_2exp(f, u, &texp, prec, rnd);

        fmpz_clear(za);
        fmpz_clear(zb);
        fmpz_clear(zc);
        fmpz_clear(zd);
        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(v);
    }
    else
    {
        ARF_MUL_TMP_DECL

        acn = an + cn;
        bdn = bn + dn;
        adn = an + dn;
        bcn = bn + cn;

        alloc = acn + bdn + adn + bcn;

        ARF_MUL_TMP_ALLOC(tmp, alloc)
        acp = tmp;
        bdp = acp + acn;
        adp = bdp + bdn;
        bcp = adp + adn;

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(acp, ap, an, cp, cn)
        acn -= (acp[0] == 0);
        acp += (acp[0] == 0);

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(bdp, bp, bn, dp, dn)
        bdn -= (bdp[0] == 0);
        bdp += (bdp[0] == 0);

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(adp, ap, an, dp, dn)
        adn -= (adp[0] == 0);
        adp += (adp[0] == 0);

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(bcp, bp, bn, cp, cn)
        bcn -= (bcp[0] == 0);
        bcp += (bcp[0] == 0);

        texp = aexp + cexp;
        uexp = bexp + dexp;
        shift = texp - uexp;

        if (shift >= 0)
            inex1 = _arf_add_mpn(e, acp, acn, asgn ^ csgn, &texp,
                                    bdp, bdn, bsgn ^ dsgn ^ 1, shift, prec, rnd);
        else
            inex1 = _arf_add_mpn(e, bdp, bdn, bsgn ^ dsgn ^ 1, &uexp,
                                    acp, acn, asgn ^ csgn, -shift, prec, rnd);

        texp = aexp + dexp;
        uexp = bexp + cexp;
        shift = texp - uexp;

        if (shift >= 0)
            inex2 = _arf_add_mpn(f, adp, adn, asgn ^ dsgn, &texp,
                                    bcp, bcn, bsgn ^ csgn, shift, prec, rnd);
        else
            inex2 = _arf_add_mpn(f, bcp, bcn, bsgn ^ csgn, &uexp,
                                    adp, adn, asgn ^ dsgn, -shift, prec, rnd);

        ARF_MUL_TMP_FREE(tmp, alloc)
    }

    return inex1 | (inex2 << 1);
}

int arf_complex_sqr(arf_t e, arf_t f,
    const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)
{
    if (!ARF_IS_LAGOM(a) || !ARF_IS_LAGOM(b) ||
        ARF_IS_SPECIAL(a) || ARF_IS_SPECIAL(b))
    {
        return arf_complex_mul_fallback(e, f, a, b, a, b, prec, rnd);
    }
    else
    {
        mp_srcptr ap, bp;
        int inex1, inex2;
        mp_ptr tmp, aap, bbp;
        mp_size_t an, bn, aan, bbn, alloc;
        slong shift;
        slong aexp, bexp;
        fmpz texp, uexp;
        TMP_INIT;

        ARF_GET_MPN_READONLY(ap, an, a);
        aexp = ARF_EXP(a);

        ARF_GET_MPN_READONLY(bp, bn, b);
        bexp = ARF_EXP(b);

        if (an >= 112 && FLINT_ABS(an - bn) <= 2 && FLINT_ABS(aexp - bexp) <= 64)
        {
            fmpzi_t x, z;
            slong abot, bbot;
            int asgn, bsgn;

            asgn = ARF_SGNBIT(a);
            bsgn = ARF_SGNBIT(b);

            abot = aexp - an * FLINT_BITS;
            bbot = bexp - bn * FLINT_BITS;

            texp = FLINT_MIN(abot, bbot);

            fmpzi_init(x);
            fmpzi_init(z);

            fmpz_lshift_mpn(fmpzi_realref(x), ap, an, asgn, abot - texp);
            fmpz_lshift_mpn(fmpzi_imagref(x), bp, bn, bsgn, bbot - texp);

            fmpzi_sqr(z, x);
            texp *= 2;

            inex1 = arf_set_round_fmpz_2exp(e, fmpzi_realref(z), &texp, prec, rnd);
            inex2 = arf_set_round_fmpz_2exp(f, fmpzi_imagref(z), &texp, prec, rnd);

            fmpzi_clear(x);
            fmpzi_clear(z);
        }
        else
        {
            aan = 2 * an;
            bbn = 2 * bn;

            alloc = aan + bbn;

            TMP_START;

            tmp = TMP_ALLOC(alloc * sizeof(mp_limb_t));
            aap = tmp;
            bbp = tmp + aan;

            FLINT_MPN_MUL_WITH_SPECIAL_CASES(aap, ap, an, ap, an)
            aan -= (aap[0] == 0);
            aap += (aap[0] == 0);

            FLINT_MPN_MUL_WITH_SPECIAL_CASES(bbp, bp, bn, bp, bn)
            bbn -= (bbp[0] == 0);
            bbp += (bbp[0] == 0);

            texp = aexp + aexp;
            uexp = bexp + bexp;
            shift = texp - uexp;

            inex2 = arf_mul(f, a, b, prec, rnd);
            ARF_EXP(f) += 1;

            if (shift >= 0)
                inex1 = _arf_add_mpn(e, aap, aan, 0, &texp,
                                        bbp, bbn, 1, shift, prec, rnd);
            else
                inex1 = _arf_add_mpn(e, bbp, bbn, 1, &uexp,
                                        aap, aan, 0, -shift, prec, rnd);

            TMP_END;
        }

        return inex1 | (inex2 << 1);
    }
}

