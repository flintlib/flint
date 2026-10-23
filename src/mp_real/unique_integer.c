/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* With the bottom limb at B^s, s = exp - size, the ball is
   M +- r with M = (d, size) B^s >= 0 up to the sign (the integers of
   -X are those of X negated) and r = err B^s.

   s >= 0: M is an integer and r >= 1 unless err = 0, so the ball holds
   a unique integer iff it is exact.

   s < 0, f = -s fractional limbs: r < 1 (err < B), so the ball holds at
   most two integers.  If size < f, then M < B^-1 and the only candidate
   is 0, contained iff M <= r, i.e. iff size = 0 or (size = 1 and
   d[0] <= err).  Otherwise M = I + F with the integer part I =
   d[f..size) and the fraction F = d[0..f) B^-f; with c the carry of
   F + err out of the f limbs (floor(M + r) = I + c), b the borrow of
   F - err (floor(M - r) = I - b) and z = [F = err] (M - r integral when
   b = 0), ceil(M - r) = I - b + (1 - z) when b = 0 (and I when b = 1,
   F - err + B^f being nonzero), and the unique integer exists iff
   floor(M + r) = ceil(M - r), i.e. iff c + b = 1 - z with b = 0, or
   c = 0 with b = 1; it is then I + c.  All three flags need only scan
   the fraction from its top until a limb decides them. */
int
mp_real_unique_integer(mp_real_t res, const mp_real_t x)
{
    slong size = x->size, s = x->exp - size, f, m, i, t;
    ulong err = x->err;
    nn_srcptr d = x->d;
    int c, b, z;

    if (s >= 0)
    {
        if (err != 0)
            return 0;
        mp_real_set(res, x);
        return 1;
    }

    f = -s;

    if (size < f)
    {
        if (!(size == 0 || (size == 1 && d[0] <= err)))
            return 0;
        mp_real_zero(res);
        return 1;
    }

    /* the flags from the limbs above the bottom one: all ones (needed for
       a carry), all zeros (needed for a borrow or equality) */
    {
        int ones = 1, zeros = 1;
        for (i = f - 1; i >= 1 && (ones || zeros); i--)
        {
            ones = ones && (d[i] == ~UWORD(0));
            zeros = zeros && (d[i] == 0);
        }
        c = ones && (d[0] + err < d[0]);
        b = zeros && (d[0] < err);
        z = zeros && (d[0] == err);
    }

    if (b ? (c != 0) : (c != 1 - z))
        return 0;

    /* n = I + c, I = (d + f, m), in exact form: its low zero limbs
       dropped, which are those of I without a carry, and those where
       I has all ones with one */
    m = size - f;
    for (t = 0; t < m && d[f + t] == (c ? ~UWORD(0) : UWORD(0)); t++)
        ;
    if (t == m)
    {
        /* n = c B^m (I = 0 without a carry, B^m - 1 with one) */
        if (c == 0)
            mp_real_zero(res);
        else
        {
            int neg = x->negative;
            mp_real_set_ui(res, 1);
            res->exp = m + 1;
            res->negative = neg;
        }
        return 1;
    }

    /* the limbs move down by f + t: in place when aliased (each limb is
       read before its new position is written), without reallocating */
    {
        int neg = x->negative;
        ulong low = d[f + t] + c;
        mp_real_fit_length(res, m - t);
        res->d[0] = low;
        flint_mpn_copyi(res->d + 1, x->d + f + t + 1, m - t - 1);
        res->size = m - t;
        res->exp = m;
        res->negative = neg;
        res->err = 0;
    }
    return 1;
}
