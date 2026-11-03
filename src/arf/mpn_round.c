/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

/*
    Rounding of limb-array mantissas (see the _arf_mpn_ interface in
    arf.h): _arf_mpn_set_round rounds an exact value, _arf_mpn_round_guarded
    an approximation with a guard limb when the guard limb certifies the
    rounding.
*/

FLINT_TLS_PREFIX nn_ptr __arf_mpn_tmp = NULL;
FLINT_TLS_PREFIX slong __arf_mpn_tmp_alloc = 0;

void
_arf_mpn_tmp_cleanup(void)
{
    flint_free(__arf_mpn_tmp);
    __arf_mpn_tmp = NULL;
    __arf_mpn_tmp_alloc = 0;
}

/* the top rn limbs of x (xn >= rn limbs) shifted left by lead bits, moved
   to res: a separate buffer, x itself, or a lower position in the same
   buffer when xn > rn (in-place rounding of a longer window) */
FLINT_FORCE_INLINE void
_move_shift(nn_ptr res, nn_srcptr x, slong xn, slong rn, unsigned int lead)
{
    nn_srcptr xs = x + xn - rn;

    if (lead == 0)
    {
        if (res < xs)
            flint_mpn_copyi(res, xs, rn);
        else if (res > xs)
            flint_mpn_copyd(res, xs, rn);
    }
    else if (res >= xs || xn == rn)
    {
        /* (mpn_lshift allows the destination at or above the source; with
           xn = rn and the destination below, the buffers are separate) */
        ulong below = (xn > rn) ? xs[-1] >> (FLINT_BITS - lead) : 0;
        FLINT_ASSERT(res >= xs || res + rn <= xs);
        mpn_lshift(res, xs, rn, lead);
        res[0] |= below;
    }
    else
    {
        /* the destination below the source, a limb of x lying below xs: a
           right shift of the window one limb lower, which mpn_rshift
           allows, the top limb's bits added */
        ulong top = xs[rn - 1] << lead;
        mpn_rshift(res, xs - 1, rn, FLINT_BITS - lead);
        res[rn - 1] |= top;
    }
}

int
_arf_mpn_set_round(nn_ptr res, slong * exp, nn_srcptr x, slong xn, slong prec, arf_rnd_t rnd)
{
    unsigned int lead;
    slong bits, rn, sb, tp, q, i;
    ulong mask;
    int inexact, increment;

    FLINT_ASSERT(prec >= 1 && xn >= 1 && x[xn - 1] != 0);

    lead = flint_clz(x[xn - 1]);
    bits = xn * FLINT_BITS - lead;
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    sb = rn * FLINT_BITS - prec;
    tp = bits - prec;               /* the discarded part is bits [0, tp) */
    *exp = -(slong) lead;

    if (tp <= 0)
    {
        /* exact: the value top-aligned in rn limbs, zero padded below */
        if (rn > xn)
        {
            _move_shift(res + rn - xn, x, xn, xn, lead);
            _flint_mpn_zero_small(res, rn - xn);
        }
        else
            _move_shift(res, x, xn, rn, lead);
        return 0;
    }

    /* the rounding decision, from x before anything is overwritten */
    q = (ulong) tp / FLINT_BITS;        /* limb of the truncation point */
    {
        slong r = (ulong) tp % FLINT_BITS;
        ulong below = (r == 0) ? 0 : (x[q] & ((UWORD(1) << r) - 1));   /* bits [q FLINT_BITS, tp) */
        int sticky;

        if (rnd == ARF_RND_NEAR)
        {
            /* the highest discarded bit, then the rest as a sticky flag */
            ulong half = (r != 0) ? (x[q] >> (r - 1)) & 1 : (x[q - 1] >> (FLINT_BITS - 1));
            ulong rest = (r != 0) ? (below & ((UWORD(1) << (r - 1)) - 1)) : (x[q - 1] & (~UWORD(0) >> 1));
            slong lim = (r != 0) ? q : q - 1;

            sticky = (rest != 0);
            for (i = lim - 1; i >= 0 && !sticky; i--)
                sticky = (x[i] != 0);

            inexact = (half != 0) || sticky;
            if (!half)
                increment = 0;
            else if (sticky)
                increment = 1;
            else    /* a tie rounds to even: the lowest kept bit is bit tp */
                increment = (x[(ulong) tp / FLINT_BITS] >> ((ulong) tp % FLINT_BITS)) & 1;
        }
        else
        {
            sticky = (below != 0);
            for (i = q - 1; i >= 0 && !sticky; i--)
                sticky = (x[i] != 0);
            inexact = sticky;
            increment = (rnd == ARF_RND_UP) && inexact;
        }
    }

    /* the kept bits, top-aligned in rn limbs */
    _move_shift(res, x, xn, rn, lead);
    mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
    res[0] &= mask;

    if (increment && _arf_mpn_increment(res, rn, sb))
        *exp += 1;

    return inexact;
}

int
_arf_mpn_round_guarded(nn_ptr res, slong * exp, nn_srcptr a, slong n, ulong g, int e, slong prec, arf_rnd_t rnd, int side)
{
    ulong mask, v;
    slong rn, sb;
    int increment;

    FLINT_ASSERT(prec >= 1 && prec <= n * FLINT_BITS && e >= 0 && e <= FLINT_BITS - 2);

    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    sb = rn * FLINT_BITS - prec;
    *exp = 0;

    /* the window: bits e .. FLINT_BITS - 2 of the guard limb, neither all
       zero nor all ones (cf. the certification argument in mul.c); with
       the value known to lie above the approximation, all zeros will do
       (the value is then just above a boundary, and rounds like the
       approximation with sticky bits), with the value below it, all ones */
    mask = (UWORD(1) << (FLINT_BITS - 1 - e)) - 1;
    v = (g >> e) & mask;

    if (v == mask && side >= 0)
        return -1;

    if (v == 0 && side <= 0)
    {
        ulong low = (sb == 0) ? 0 : (a[n - rn] & ((UWORD(1) << sb) - 1));

        if (side == 0 || g != 0 || low != 0 || !flint_mpn_zero_p(a, n - rn))
            return -1;

        /* the value lies just below the approximation, which is exact at
           prec bits: the rounding is decided (it is never a tie) */
        if (res != a + n - rn)
            flint_mpn_copyi(res, a + n - rn, rn);
        if (rnd == ARF_RND_DOWN)
            *exp = _arf_mpn_decrement(res, rn, sb);
        return 1;
    }

    if (rnd == ARF_RND_DOWN)
        increment = 0;
    else if (rnd == ARF_RND_UP)
        increment = 1;
    else
    {
        /* nearest: the highest discarded bit (no tie is possible) */
        if (sb != 0)
            increment = (a[n - rn] >> (sb - 1)) & 1;
        else if (rn < n)
            increment = a[n - rn - 1] >> (FLINT_BITS - 1);
        else
            increment = g >> (FLINT_BITS - 1);
    }

    if (res != a + n - rn)
        flint_mpn_copyi(res, a + n - rn, rn);
    mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
    res[0] &= mask;

    if (increment && _arf_mpn_increment(res, rn, sb))
        *exp = 1;

    return 1;
}

/* Is rounding {a, an} (an approximation within 2^e units of a[0] of the
   value P, on either side) to prec bits guaranteed to give the same
   result as rounding P, in every strict rounding mode, with P inexact?
   The bits of a from bit e of its lowest limb up to the bit below the
   rounding point (at most one limb of them) must be neither all zero
   nor all ones. */
int
_arf_rounding_certain(nn_srcptr a, slong an, slong e, slong prec)
{
    slong k, w, q, r;
    ulong v, mask;

    /* number of bits of A below the rounding point */
    k = an * FLINT_BITS - flint_clz(a[an - 1]) - prec;
    /* bits [e, k-1) are available; inspect at most one limb of them */
    w = k - 1 - e;

    if (w < 4)
        return 0;

    w = FLINT_MIN(w, FLINT_BITS);

    q = (ulong) e / FLINT_BITS;
    r = (ulong) e % FLINT_BITS;

    v = a[q] >> r;
    if (r != 0 && r + w > FLINT_BITS)
        v |= a[q + 1] << (FLINT_BITS - r);

    mask = (w == FLINT_BITS) ? ~UWORD(0) : ((UWORD(1) << w) - 1);
    v &= mask;

    return (v != 0) && (v != mask);
}
