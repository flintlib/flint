/*
    Copyright (C) 2014, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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



/* n limbs of |x| 2^-shift (0 <= shift < FLINT_BITS) aligned to the top of
   an n-limb frame, truncating below and zero-padding; scratch has n + 2
   limbs */
static void
_arf_complex_part(nn_ptr dst, nn_srcptr xp, slong xn, slong n, slong shift, nn_ptr scratch)
{
    slong m = FLINT_MIN(xn, n + 1);

    if (shift == 0)
    {
        if (xn >= n)
        {
            flint_mpn_copyi(dst, xp + xn - n, n);
        }
        else
        {
            flint_mpn_zero(dst, n - xn);
            flint_mpn_copyi(dst + n - xn, xp, xn);
        }
        return;
    }

    flint_mpn_zero(scratch, n + 1 - m);
    flint_mpn_copyi(scratch + n + 1 - m, xp + xn - m, m);
    mpn_rshift(scratch, scratch, n + 1, shift);
    flint_mpn_copyi(dst, scratch + 1, n);
}

/*
    Approximate complex product e + f i = (a + b i)(c + d i) for use in dot
    products, where only the bits above 2^bottom matter: the four parts
    are at least n limbs long (n as below) and balanced,
    the parts are truncated to n limbs at the top exponent of each operand,
    flint_mpn_mulhigh_n_complex returns the limbs [n, 2n] of each product
    with an error below 3 ulp of its lowest limb, and the four truncations
    add less than 4 ulp. With n chosen so that this lowest limb lies at
    least 64 bits below the bottom of both accumulators, the total error
    is below one accumulator ulp for each output, which the caller charges
    to the error counters. e and f receive the (unrounded) approximations.
    Returns 0 if the shape is unsuitable (the parts are not all at least
    3n/4 limbs long, or the exact products would barely exceed n limbs).
*/
int
_arf_complex_mul_high_approx(arf_t e, arf_t f, const arf_t a, const arf_t b,
    const arf_t c, const arf_t d, slong bottom)
{
    nn_srcptr ap, bp, cp, dp;
    slong an, bn, cn, dn, n, E1, E2, term_prec, alloc, sn, fix, i;
    nn_ptr tmp, A, B, C, D, zr, zi, scratch;
    int zr_sgn, zi_sgn;
    TMP_INIT;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);

    E1 = FLINT_MAX(ARF_EXP(a), ARF_EXP(b));
    E2 = FLINT_MAX(ARF_EXP(c), ARF_EXP(d));

    /* bits of the products above the bottom, plus one limb */
    term_prec = (E1 + E2) - bottom + FLINT_BITS;
    if (term_prec <= 0)
        return 0;
    n = (term_prec + FLINT_BITS - 1) / FLINT_BITS;

    /* all parts at least n limbs (no padding: in a dot product the
       selection must stay cheap for the many terms it does not apply to),
       and the kernel size from which it beats four real products */
    if (n < 13 || n > 2000 ||
        an < n || bn < n || cn < n || dn < n ||
        FLINT_ABS(ARF_EXP(a) - ARF_EXP(b)) >= FLINT_BITS ||
        FLINT_ABS(ARF_EXP(c) - ARF_EXP(d)) >= FLINT_BITS)
        return 0;

    alloc = 4 * n + 2 * (n + 1) + (n + 2);
    TMP_START;
    tmp = TMP_ALLOC(sizeof(ulong) * alloc);
    A = tmp; B = A + n; C = B + n; D = C + n;
    zr = D + n; zi = zr + n + 1; scratch = zi + n + 1;

    _arf_complex_part(A, ap, an, n, E1 - ARF_EXP(a), scratch);
    _arf_complex_part(B, bp, bn, n, E1 - ARF_EXP(b), scratch);
    _arf_complex_part(C, cp, cn, n, E2 - ARF_EXP(c), scratch);
    _arf_complex_part(D, dp, dn, n, E2 - ARF_EXP(d), scratch);

    flint_mpn_mulhigh_n_complex(zr, &zr_sgn, zi, &zi_sgn,
        A, ARF_SGNBIT(a), B, ARF_SGNBIT(b), C, ARF_SGNBIT(c), D, ARF_SGNBIT(d), n);

    for (i = 0; i < 2; i++)
    {
        nn_ptr z = (i == 0) ? zr : zi;
        int sgn = (i == 0) ? zr_sgn : zi_sgn;
        arf_ptr out = (i == 0) ? e : f;

        sn = n + 1;
        while (sn > 0 && z[sn - 1] == 0)
            sn--;

        if (sn == 0)
        {
            arf_zero(out);
        }
        else
        {
            /* value = z / B^(n+1) * 2^(E1 + E2 + FLINT_BITS) */
            _arf_set_round_mpn(out, &fix, z, sn, sgn, ARF_PREC_EXACT, ARF_RND_DOWN);
            fmpz_set_si(ARF_EXPREF(out), E1 + E2 + FLINT_BITS + fix + (sn - (n + 1)) * FLINT_BITS);
        }
    }

    TMP_END;
    return 1;
}



/* dst = x << shift as an integer, returning its length (at most
   xn + shift / FLINT_BITS + 1 limbs, possibly with a zero top limb) */

static slong
_arf_lshift_limbs(nn_ptr dst, nn_srcptr xp, slong xn, slong shift)
{
    slong limbs = shift / FLINT_BITS, bits = shift % FLINT_BITS;

    flint_mpn_zero(dst, limbs);

    if (bits == 0)
    {
        flint_mpn_copyi(dst + limbs, xp, xn);
        return xn + limbs;
    }
    else
    {
        dst[limbs + xn] = mpn_lshift(dst + limbs, xp, xn, bits);
        return xn + limbs + 1;
    }
}

/*
    Exact complex product e + f i = (a + b i)(c + d i) of finite, nonzero,
    LAGOM operands via flint_mpn_mul_complex: the two parts of each
    operand are aligned to a common bottom exponent (the shorter
    exponent gap must be small for this to be efficient; the caller
    checks the shape).
*/
void
_arf_complex_mul_exact(arf_t e, arf_t f, const arf_t a, const arf_t b,
    const arf_t c, const arf_t d)
{
    nn_srcptr ap, bp, cp, dp;
    slong an, bn, cn, dn, An, Bn, Cn, Dn, abot, bbot, cbot, dbot, texp, uexp;
    slong alloc, outn, lr, li, fix;
    nn_ptr tmp, A, B, C, D, zr, zi;
    int sqr;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);

    abot = ARF_EXP(a) - an * FLINT_BITS;
    bbot = ARF_EXP(b) - bn * FLINT_BITS;
    cbot = ARF_EXP(c) - cn * FLINT_BITS;
    dbot = ARF_EXP(d) - dn * FLINT_BITS;
    texp = FLINT_MIN(abot, bbot);
    uexp = FLINT_MIN(cbot, dbot);

    sqr = (a == c && b == d);

    An = an + (abot - texp) / FLINT_BITS + 1;
    Bn = bn + (bbot - texp) / FLINT_BITS + 1;
    Cn = cn + (cbot - uexp) / FLINT_BITS + 1;
    Dn = dn + (dbot - uexp) / FLINT_BITS + 1;
    outn = FLINT_MAX(An, Bn) + FLINT_MAX(Cn, Dn) + 1;
    alloc = An + Bn + Cn + Dn + 2 * outn;
    tmp = flint_malloc(sizeof(ulong) * alloc);
    A = tmp; B = A + An; C = B + Bn; D = C + Cn; zr = D + Dn; zi = zr + outn;

    An = _arf_lshift_limbs(A, ap, an, abot - texp);
    Bn = _arf_lshift_limbs(B, bp, bn, bbot - texp);

    if (sqr)
    {
        flint_mpn_sqr_complex(zr, &lr, zi, &li, A, An, ARF_SGNBIT(a), B, Bn, ARF_SGNBIT(b));
    }
    else
    {
        Cn = _arf_lshift_limbs(C, cp, cn, cbot - uexp);
        Dn = _arf_lshift_limbs(D, dp, dn, dbot - uexp);
        flint_mpn_mul_complex(zr, &lr, zi, &li, A, An, ARF_SGNBIT(a), B, Bn, ARF_SGNBIT(b),
            C, Cn, ARF_SGNBIT(c), D, Dn, ARF_SGNBIT(d));
    }

    /* value = z 2^(texp + uexp) */
    if (lr == 0)
        arf_zero(e);
    else
    {
        slong n = FLINT_ABS(lr);
        while (zr[n - 1] == 0) n--;
        _arf_set_round_mpn(e, &fix, zr, n, lr < 0, ARF_PREC_EXACT, ARF_RND_DOWN);
        fmpz_set_si(ARF_EXPREF(e), texp + uexp + n * FLINT_BITS + fix);
    }

    if (li == 0)
        arf_zero(f);
    else
    {
        slong n = FLINT_ABS(li);
        while (zi[n - 1] == 0) n--;
        _arf_set_round_mpn(f, &fix, zi, n, li < 0, ARF_PREC_EXACT, ARF_RND_DOWN);
        fmpz_set_si(ARF_EXPREF(f), texp + uexp + n * FLINT_BITS + fix);
    }

    flint_free(tmp);
}

/* precisions for which _arf_complex_mul_high can apply: at least 13
   output limbs (below, four real products on the stack are as fast as
   the complex high product kernel; measured) and not ARF_PREC_EXACT */
#define ARF_COMPLEX_MUL_HIGH_PREC_OK(prec) \
    ((prec) > 12 * FLINT_BITS && (prec) != ARF_PREC_EXACT)

static int
_arf_complex_mul_high(arf_t e, arf_t f, const arf_t a, const arf_t b,
    const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)
{
    nn_srcptr ap, bp, cp, dp;
    slong an, bn, cn, dn, n, outn, E1, E2, alloc, sn, fix, k, i;
    nn_ptr tmp, A, B, C, D, zr, zi, scratch;
    int zr_sgn, zi_sgn, relaxed, inex[2], sqr;
    arf_ptr out[2];
    TMP_INIT;

    outn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    n = outn + (FLINT_BITS >= 64 ? 1 : 2);

    /* balanced parts of comparable magnitude, and long enough that the
       kernel beats four real products (which also truncate); tested on
       the size fields first, as this fails for most calls at moderate
       precision (the precision itself is checked by the caller, before
       the call) */
    FLINT_ASSERT(ARF_COMPLEX_MUL_HIGH_PREC_OK(prec));
    if (n > 2000 ||
        4 * ARF_SIZE(a) < 3 * n || 4 * ARF_SIZE(b) < 3 * n ||
        4 * ARF_SIZE(c) < 3 * n || 4 * ARF_SIZE(d) < 3 * n ||
        FLINT_ABS(ARF_EXP(a) - ARF_EXP(b)) >= FLINT_BITS ||
        FLINT_ABS(ARF_EXP(c) - ARF_EXP(d)) >= FLINT_BITS)
        return -1;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);

    /* the exact products are short enough to prefer */
    if (an + cn <= outn + FLINT_MAX(2, outn / 16))
        return -1;

    relaxed = (rnd >= ARF_RND_FAST);
    sqr = (a == c && b == d);

    E1 = FLINT_MAX(ARF_EXP(a), ARF_EXP(b));
    E2 = FLINT_MAX(ARF_EXP(c), ARF_EXP(d));

    alloc = 4 * n + 2 * (n + 1) + (n + 2);
    TMP_START;
    tmp = TMP_ALLOC(sizeof(ulong) * alloc);
    A = tmp; B = A + n; C = B + n; D = C + n;
    zr = D + n; zi = zr + n + 1; scratch = zi + n + 1;

    _arf_complex_part(A, ap, an, n, E1 - ARF_EXP(a), scratch);
    _arf_complex_part(B, bp, bn, n, E1 - ARF_EXP(b), scratch);

    if (sqr)
    {
        flint_mpn_sqrhigh_n_complex(zr, &zr_sgn, zi, &zi_sgn,
            A, ARF_SGNBIT(a), B, ARF_SGNBIT(b), n);
    }
    else
    {
        _arf_complex_part(C, cp, cn, n, E2 - ARF_EXP(c), scratch);
        _arf_complex_part(D, dp, dn, n, E2 - ARF_EXP(d), scratch);
        flint_mpn_mulhigh_n_complex(zr, &zr_sgn, zi, &zi_sgn,
            A, ARF_SGNBIT(a), B, ARF_SGNBIT(b), C, ARF_SGNBIT(c), D, ARF_SGNBIT(d), n);
    }

    out[0] = e; out[1] = f;
    inex[0] = inex[1] = -1;

    for (i = 0; i < 2; i++)
    {
        nn_ptr z = (i == 0) ? zr : zi;

        sn = n + 1;
        while (sn > 0 && z[sn - 1] == 0)
            sn--;

        if (sn == 0)
            break;

        k = sn * FLINT_BITS - flint_clz(z[sn - 1]) - prec;

        if (relaxed)
        {
            if (k - 1 - 4 < 30)
                break;
            /* both outputs must succeed: defer writing until checked */
            inex[i] = 1;
        }
        else
        {
            if (!_arf_rounding_certain(z, sn, 4, prec))
                break;
            inex[i] = 0;
        }
    }

    if (inex[0] == -1 || inex[1] == -1)
    {
        TMP_END;
        return -1;
    }

    for (i = 0; i < 2; i++)
    {
        nn_ptr z = (i == 0) ? zr : zi;
        int sgn = (i == 0) ? zr_sgn : zi_sgn;

        sn = n + 1;
        while (z[sn - 1] == 0)
            sn--;

        inex[i] = _arf_set_round_mpn(out[i], &fix, z, sn, sgn, prec,
            relaxed ? ARF_RND_NEAR : rnd);
        if (relaxed)
            inex[i] = 1;
        /* value = z / B^(n+1) * 2^(E1 + E2 + FLINT_BITS) */
        fmpz_set_si(ARF_EXPREF(out[i]), E1 + E2 + FLINT_BITS + fix + (sn - (n + 1)) * FLINT_BITS);
    }

    TMP_END;
    return inex[0] | (inex[1] << 1);
}

int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b,
                                      const arf_t c, const arf_t d,
                                      slong prec, arf_rnd_t rnd)
{
    nn_srcptr ap, bp, cp, dp;
    int asgn, bsgn, csgn, dsgn, inex1, inex2;
    nn_ptr tmp, acp, bdp, adp, bcp;
    slong an, bn, cn, dn, acn, bdn, adn, bcn, alloc;
    slong shift;
    slong aexp, bexp, cexp, dexp;
    fmpz texp, uexp;

    if (!ARF_IS_LAGOM(a) || !ARF_IS_LAGOM(b) ||
        !ARF_IS_LAGOM(c) || !ARF_IS_LAGOM(d) ||
        ARF_IS_SPECIAL(d) || ARF_IS_SPECIAL(a) ||
        ARF_IS_SPECIAL(b) || ARF_IS_SPECIAL(c))
    {
        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);
        return arf_complex_mul_fallback(e, f, a, b, c, d, prec, rnd);
    }

    if (ARF_COMPLEX_MUL_HIGH_PREC_OK(prec))
    {
        int r = _arf_complex_mul_high(e, f, a, b, c, d, prec, rnd);
        if (r != -1)
            return r;
    }

    /* the remaining paths use strict rounding only */
    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);

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

    /* exact complex product (Karatsuba / FFT inside the kernel) */
    if (an >= 20 &&
        cn >= 20 &&
        FLINT_ABS(an - bn) <= 2 &&
        FLINT_ABS(cn - dn) <= 2 &&
        FLINT_ABS(aexp - bexp) <= 64 &&
        FLINT_ABS(cexp - dexp) <= 64)
    {
        _arf_complex_mul_exact(e, f, a, b, c, d);
        inex1 = arf_set_round(e, e, prec, rnd);
        inex2 = arf_set_round(f, f, prec, rnd);
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

/* the exact complex square (flint_mpn_sqr_complex) from this many limbs,
   below two real squares and a product (measured on Skylake: 0.95 of the
   time at 25 limbs, 0.8-0.9 at 30-200, 0.7 from 500) */
#define ARF_COMPLEX_SQR_EXACT_CUTOFF 24

int arf_complex_sqr(arf_t e, arf_t f,
    const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)
{
    if (!ARF_IS_LAGOM(a) || !ARF_IS_LAGOM(b) ||
        ARF_IS_SPECIAL(a) || ARF_IS_SPECIAL(b))
    {
        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);
        return arf_complex_mul_fallback(e, f, a, b, a, b, prec, rnd);
    }

    if (ARF_COMPLEX_MUL_HIGH_PREC_OK(prec))
    {
        int r = _arf_complex_mul_high(e, f, a, b, a, b, prec, rnd);
        if (r != -1)
            return r;
    }

    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);

    {
        nn_srcptr ap, bp;
        int inex1, inex2;
        nn_ptr tmp, aap, bbp;
        slong an, bn, aan, bbn, alloc;
        slong shift;
        slong aexp, bexp;
        fmpz texp, uexp;
        TMP_INIT;

        ARF_GET_MPN_READONLY(ap, an, a);
        aexp = ARF_EXP(a);

        ARF_GET_MPN_READONLY(bp, bn, b);
        bexp = ARF_EXP(b);

        if (an >= ARF_COMPLEX_SQR_EXACT_CUTOFF && FLINT_ABS(an - bn) <= 2 && FLINT_ABS(aexp - bexp) <= 64)
        {
            /* the exact complex square (flint_mpn_sqr_complex) */
            _arf_complex_mul_exact(e, f, a, b, a, b);
            inex1 = arf_set_round(e, e, prec, rnd);
            inex2 = arf_set_round(f, f, prec, rnd);
        }
        else
        {
            aan = 2 * an;
            bbn = 2 * bn;

            alloc = aan + bbn;

            TMP_START;

            tmp = TMP_ALLOC(alloc * sizeof(ulong));
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
