/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* The slow paths of the radius normalization (see impl.h). */

/* strip high zero limbs (and low zero limbs when exact); zero
   mantissas collapse to size 0 keeping the anchor in exp */
void
_mp_real_norm_slow(mp_real_t x)
{
    while (x->size > 0 && x->d[x->size - 1] == 0)
    {
        x->size--;
        x->exp--;
    }

    if (x->size == 0)
    {
        x->negative = 0;
        return;
    }

    if (x->err == 0)
    {
        slong t;
        for (t = 0; x->d[t] == 0; t++)
            ;
        if (t > 0)
        {
            flint_mpn_copyi(x->d, x->d + t, x->size - t);
            x->size -= t;
        }
    }
}

/* Install the composed bound e on x (whose mantissa, sign and
   exp/size are already set and top-normalized), as a count of ulps of
   the bottom limb: a count of B or more truncates the mantissa by as
   many limbs as needed (each dropped run of limbs is worth one unit),
   a bound below one ulp pads the mantissa with zero limbs down to its
   scale. */
void
_mp_real_apply_bnd_slow(mp_real_t x, mp_real_bnd_t e)
{
    slong anc, k;
    ulong v;

    v = _mp_real_bnd_reduce(&e);

    if (x->size == 0)
    {
        /* a zero midpoint anchors its radius at exp: no padding */
        x->exp = e.a;
        x->err = v;
        return;
    }

    anc = x->exp - x->size;
    k = e.a - anc;          /* e = v B^k units of the bottom limb */

    if (k > 0)
    {
        if (k > x->size - 1)
        {
            /* the bound swamps the mantissa: 0 +/- (|x| + e) with
               |x| < B^exp, at the anchor B^exp or above */
            mp_real_bnd_t t = _mp_real_bnd_add(_mp_real_bnd(0, v, e.a), _mp_real_bnd(0, 1, x->exp));
            v = _mp_real_bnd_reduce(&t);
            x->size = 0;
            x->negative = 0;
            x->exp = t.a;
            x->err = v;
            return;
        }
        /* drop k limbs, worth one unit */
        flint_mpn_copyi(x->d, x->d + k, x->size - k);
        x->size -= k;
        v++;
        if (v == 0)
        {
            /* the count was B - 1: one more limb, two units */
            flint_mpn_copyi(x->d, x->d + 1, x->size - 1);
            x->size -= 1;
            v = 2;
        }
    }
    else if (k < 0)
    {
        /* pad -k zero limbs below the mantissa */
        mp_real_fit_length(x, x->size - k);
        memmove(x->d - k, x->d, x->size * sizeof(ulong));
        flint_mpn_zero(x->d, -k);
        x->size -= k;
    }

    x->err = v;

    if (x->size > 0 && x->d[x->size - 1] == 0)
        _mp_real_norm(x);
}
