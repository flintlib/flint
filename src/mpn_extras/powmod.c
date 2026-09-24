/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#define EBIT(e, i) (((e)[(i) / FLINT_BITS] >> ((i) % FLINT_BITS)) & 1)

/* Window width, cut back so that the table of odd powers stays modest. */
static int
_window_size(flint_bitcnt_t ebits, mp_size_t n)
{
    int w;

    if (ebits < 8)
        w = 1;
    else if (ebits < 24)
        w = 2;
    else if (ebits < 70)
        w = 3;
    else if (ebits < 200)
        w = 4;
    else if (ebits < 500)
        w = 5;
    else if (ebits < 1200)
        w = 6;
    else if (ebits < 3000)
        w = 7;
    else
        w = 8;

    while (w > 1 && (((mp_size_t) 1) << (w - 1)) * n > 65536)
        w--;

    return w;
}

void
flint_mpn_powmod_preinvn(mp_ptr res, mp_srcptr a, mp_srcptr e, mp_size_t en,
        mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    flint_bitcnt_t ebits;
    slong i, l, k, tabn;
    int w, started = 0;
    ulong val;
    mp_ptr tab, sqr;
    TMP_INIT;

    while (en > 0 && e[en - 1] == 0)
        en--;

    if (en == 0)
    {
        flint_mpn_zero(res, n);
        res[0] = UWORD(1) << norm;
        return;
    }

    ebits = en * FLINT_BITS - flint_clz(e[en - 1]);
    w = _window_size(ebits, n);
    tabn = ((mp_size_t) 1) << (w - 1);

    TMP_START;

    tab = TMP_ALLOC((tabn + 1) * n * sizeof(mp_limb_t));
    sqr = tab + tabn * n;

    /* the odd powers a, a^3, ..., a^(2^w - 1) */
    flint_mpn_copyi(tab, a, n);

    if (tabn > 1)
    {
        flint_mpn_mulmod_preinvn(sqr, a, a, n, d, dinv, norm);

        for (k = 1; k < tabn; k++)
            flint_mpn_mulmod_preinvn(tab + k * n, tab + (k - 1) * n, sqr,
                    n, d, dinv, norm);
    }

    for (i = ebits - 1; i >= 0; )
    {
        if (!EBIT(e, i))
        {
            if (started)
                flint_mpn_mulmod_preinvn(res, res, res, n, d, dinv, norm);

            i--;
            continue;
        }

        /* the longest window ending in a set bit, hence of odd value */
        l = FLINT_MAX(i - w + 1, 0);

        while (!EBIT(e, l))
            l++;

        val = 0;
        for (k = i; k >= l; k--)
            val = 2 * val + EBIT(e, k);

        if (started)
        {
            for (k = 0; k <= i - l; k++)
                flint_mpn_mulmod_preinvn(res, res, res, n, d, dinv, norm);

            flint_mpn_mulmod_preinvn(res, res, tab + ((val - 1) / 2) * n,
                    n, d, dinv, norm);
        }
        else
        {
            /* the leading window, which needs no squarings before it */
            flint_mpn_copyi(res, tab + ((val - 1) / 2) * n, n);
            started = 1;
        }

        i = l - 1;
    }

    TMP_END;
}
