/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "impl.h"

/*
    Hyperbolic functions of floating-point mantissas (mp_real.h) from
    E = exp |x| and its reciprocal: cosh = (E + 1/E) / 2, sinh = (E - 1/E) / 2,
    tanh = (E - 1/E) / (E + 1/E), all in fixed point.  For |x| < 1 the
    difference E - 1/E cancels about -e + 1 bits, so the exponential runs
    at as many more limbs; one approximate reciprocal (flint_mpn_invapprox)
    gives 1/E.
*/

/* Y = (y, w + 1) and R = (r, w + 1) with exp |x| = Y 2^(k - FLINT_BITS w)
   and exp(-|x|) = R 2^(k - FLINT_BITS w), both within the returned
   number of units; k = 0 for |x| < 1 (then Y < e B^w), else k >= 1 with
   Y within its error of [B^w, 2 B^w).  The reciprocal: R0 = floor(B^(2w) / Y) or one
   more, within ey + 3 units of B^(2w) / E' for the true E' = exp |x|
   2^(FLINT_BITS w - k) (|B^(2w) / Y - B^(2w) / E'| = B^(2w) |E' - Y| /
   (Y E') < ey + 1 for Y, E' >= B^w - ey), then shifted right by 2k bits
   (two more units). */
static slong
_exp_pm(nn_ptr y, nn_ptr r, ulong * err, nn_srcptr d, slong n, slong e, slong w)
{
    nn_ptr R0;
    slong k, yn;
    ulong ey, er;
    TMP_INIT;

    k = _mp_real_exp_fp(y, &ey, d, n, e, 0, FLINT_BITS * w);
    FLINT_ASSERT(k >= 0);

    TMP_START;
    R0 = TMP_ALLOC(sizeof(ulong) * (w + 2));

    yn = w + 1;
    while (y[yn - 1] == 0)
        yn--;
    /* (R0, 2 w - yn + 2), the top limb zero when yn = w */
    flint_mpn_invapprox(R0, y, yn, 2 * w);
    er = ey + 3;

    if (k == 0)
        flint_mpn_copyi(r, R0, w + 1);
    else
    {
        slong ls = (2 * k) / FLINT_BITS;
        unsigned int bs = (2 * k) % FLINT_BITS;

        if (ls >= w + 1)
            _flint_mpn_zero_small(r, w + 1);
        else
        {
            if (bs != 0)
                mpn_rshift(r, R0 + ls, w + 1 - ls, bs);
            else
                flint_mpn_copyi(r, R0 + ls, w + 1 - ls);
            _flint_mpn_zero_small(r + w + 1 - ls, ls);
        }
        er = (2 * k >= FLINT_BITS) ? 2 : (er >> (2 * k)) + 2;
    }

    *err = FLINT_MAX(ey, er);
    TMP_END;
    return k;
}

/* the limbs for |x| < 1: the -e + 1 bits cancelling in E - 1/E */
static slong
_limbs(slong pb, slong e)
{
    return (e <= 0) ? _mp_real_elem_limbs(pb + 1 - e) : _mp_real_elem_limbs(pb);
}

slong
_mp_real_sinh_cosh_fp(nn_ptr ys, nn_ptr yc, slong * wout, ulong * err, nn_srcptr d, slong n, slong e, slong pb)
{
    nn_ptr y, r;
    slong w, k;
    ulong ee;
    TMP_INIT;

    w = _limbs(pb, e);
    *wout = w;

    TMP_START;
    y = TMP_ALLOC(sizeof(ulong) * (2 * w + 4));
    r = y + w + 2;

    k = _exp_pm(y, r, &ee, d, n, e, w);

    /* (E -+ 1/E) / 2: the unit halved */
    if (ys != NULL)
        mpn_sub_n(ys, y, r, w + 1);
    if (yc != NULL)
        mpn_add_n(yc, y, r, w + 1);

    *err = 2 * ee;
    TMP_END;
    return k - 1;
}

void
_mp_real_tanh_fp(nn_ptr q, slong * wout, ulong * err, nn_srcptr d, slong n, slong e, slong pb)
{
    nn_ptr y, r, t, s;
    slong w, k;
    ulong ee;
    TMP_INIT;

    w = _limbs(pb, e);
    *wout = w;

    /* large x: tanh |x| = 1 - 2/(E^2 + 1) > 1 - 2 e^(-2|x|), one within
       a unit once 2 |x| log2(e) - 1 >= FLINT_BITS w + 1 (|x| from the top
       limb in double arithmetic, rounded down) */
    if (e >= 1 && (e > 62 || 2.8853900817779268 * (1 - 1e-15) * ldexp((double) d[n - 1], (int) (e - FLINT_BITS)) >= (double) (FLINT_BITS * w + 2)))
    {
        _flint_mpn_zero_small(q, w);
        q[w] = 1;
        *err = 1;
        return;
    }

    TMP_START;
    y = TMP_ALLOC(sizeof(ulong) * (4 * w + 8));
    r = y + w + 2;
    t = r + w + 2;
    s = t + w + 2;

    k = _exp_pm(y, r, &ee, d, n, e, w);
    (void) k;

    /* T = Y - R < S = Y + R with S >= B^w (Y >= B^w for k >= 1, Y + R >=
       2 B^w for k = 0): Q = floor(T B^w / S) or one more, and as
       d(T / S) = 2 (R dY - Y dR) / S^2, the errors ee of Y and R move
       T B^w / S by at most 2 ee B^w (Y + R) / S^2 <= 2 ee, plus two */
    mpn_sub_n(t, y, r, w + 1);
    mpn_add_n(s, y, r, w + 1);
    flint_mpn_divapprox_fraction(q, t, w + 1, s, w + 1, w);
    *err = 2 * ee + 2;
    TMP_END;
}
