/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* exp of a fraction 0 <= x < 2^-32, (res, n + 1) with a units limb.

   Dispatch: hardcoded straight-line rectangular splitting on 64-bit
   machines ("Fast multiple precision exp(x) with precomputations", van
   der Hoeven and Johansson; "Efficient implementation of elementary
   functions in the medium-precision range", Johansson) for
   2^-64 <= x < 2^-32 with n <= 11 (the table is indexed by the number
   of Taylor terms N, n = ceil(N/2); N = 21 serves n = 11 since
   21! > 2^33 damps the half-limb overhang) and for
   2^-128 <= x < 2^-64 with n <= 21; otherwise the tapered rectangular
   splitting of series_rs.c at the argument's actual leading zero bits
   (1.3 to 1.8 times faster than the untapered generic routine it
   replaced from 14 limbs up).

   Error bound (10 ulps returned): hardcoded routines <= 8 ulps plus a
   sub-ulp tail, series_rs.c within 5 ulps; all truncations are
   downward, so the bound is one-sided. */

#if FLINT_BITS == 64
#include "exp_rs_hard.inc"
#endif

void
_mp_real_exp_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    slong zb;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

    if (err != NULL)
        *err = 10;

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 11)
        {
            _mp_real_exp_rs32_tab[FLINT_MIN(2 * n, 21)](res, x);
            return;
        }
    }
    else if (n <= 21 && !(n >= 2 && x[n - 2] == 0))
    {
        /* 2^-128 <= x < 2^-64 */
        _mp_real_exp_rs_tab[n](res, x);
        return;
    }
#endif

    zb = _mp_real_elem_lzb(x, n);
    if (zb == WORD_MAX)
    {
        flint_mpn_zero(res, n);
        res[n] = 1;
    }
    else
        _mp_real_series_rs(res, x, n, zb, MP_REAL_SERIES_EXP);
}
