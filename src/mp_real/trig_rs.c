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

/* sin, cos, sinh, cosh, atan and atanh of a fraction 0 <= x < 2^-32,
   (res, n + 1) with a units limb (res, n for atan and atanh).

   Dispatch: hardcoded straight-line rectangular splitting on 64-bit
   machines for 2^-64 <= x < 2^-32 (n <= 10 for the sin/cos families,
   n <= 17 for atan/atanh) and for 2^-128 <= x < 2^-64 (n <= 22 resp.
   n <= 35), where they are fastest; otherwise the tapered rectangular
   splitting of series_rs.c at the argument's actual leading zero bits
   (1.3 to 1.7 times faster than the untapered generic routines it
   replaced, which tapered only by whole zero limbs).

   Error bounds (15 ulps returned): hardcoded routines <= 12 ulps;
   series_rs.c within 2 ulps, plus 1 for forming cos = 1 - (1 - cos).
   Alternating-series bounds are two-sided, hyperbolic ones
   one-sided. */

#if FLINT_BITS == 64
#include "trig_rs_hard.inc"
#include "hand_mulhi.inc"
#endif

/* ysin = sin x (resp. sinh x), ycos = cos x (resp. cosh x), each n + 1
   limbs, by series_rs.c */
static void
_mp_real_sc_rs_tapered(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n,
    int alternating)
{
    slong zb = _mp_real_elem_lzb(x, n);

    if (zb == WORD_MAX)
    {
        if (ysin != NULL)
            flint_mpn_zero(ysin, n + 1);
        if (ycos != NULL)
        {
            flint_mpn_zero(ycos, n);
            ycos[n] = 1;
        }
        return;
    }

    _mp_real_series_rs_sin_cos(ysin, ycos, x, n, zb, !alternating);

    if (ysin != NULL)
        ysin[n] = 0;

    if (ycos != NULL)
    {
        /* ycos holds g = 1 - cos x resp. cosh x - 1 */
        if (!alternating)
            ycos[n] = 1;
        else if (flint_mpn_zero_p(ycos, n))
            ycos[n] = 1;
        else
        {
            mpn_neg(ycos, ycos, n);
            ycos[n] = 0;
        }
    }
}

static void
_mp_real_at_rs_tapered(nn_ptr res, nn_srcptr x, slong n, int alternating)
{
    slong zb = _mp_real_elem_lzb(x, n);

    if (zb == WORD_MAX)
        flint_mpn_zero(res, n);
    else
        _mp_real_series_rs(res, x, n, zb, alternating
            ? MP_REAL_SERIES_ATAN : MP_REAL_SERIES_ATANH);
}

/* dispatch helpers ********************************************************/

static void
_mp_real_sc_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n,
    int alternating)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 10)
        {
            (alternating ? _mp_real_sin_cos_rs32_tab
                         : _mp_real_sinh_cosh_rs32_tab)[n](ysin, ycos, x);
            return;
        }
    }
    else
    {
        if (n <= 22 && !(n >= 2 && x[n - 2] == 0))
        {
            /* 2^-128 <= x < 2^-64 */
            (alternating ? _mp_real_sin_cos_rs_tab
                         : _mp_real_sinh_cosh_rs_tab)[n](ysin, ycos, x);
            return;
        }
    }
#endif
    _mp_real_sc_rs_tapered(ysin, ycos, x, n, alternating);
}

static void
_mp_real_at_rs(nn_ptr res, nn_srcptr x, slong n, int alternating)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT((x[n - 1] >> (FLINT_BITS - 32)) == 0);

#if FLINT_BITS == 64
    if (x[n - 1] != 0)
    {
        /* 2^-64 <= x < 2^-32 */
        if (n <= 17)
        {
            (alternating ? _mp_real_atan_rs32_tab
                         : _mp_real_atanh_rs32_tab)[n](res, x);
            return;
        }
    }
    else
    {
        if (n <= 35 && !(n >= 2 && x[n - 2] == 0))
        {
            /* 2^-128 <= x < 2^-64 */
            (alternating ? _mp_real_atan_rs_tab
                         : _mp_real_atanh_rs_tab)[n](res, x);
            return;
        }
    }
#endif
    _mp_real_at_rs_tapered(res, x, n, alternating);
}

/* public functions ********************************************************/

void _mp_real_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(ysin, ycos, x, n, 1);
    if (err != NULL)
        *err = 15;
}

void _mp_real_sin_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(res, NULL, x, n, 1);
    if (err != NULL)
        *err = 15;
}

void _mp_real_cos_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(NULL, res, x, n, 1);
    if (err != NULL)
        *err = 15;
}

void _mp_real_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(ysinh, ycosh, x, n, 0);
    if (err != NULL)
        *err = 15;
}

void _mp_real_sinh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(res, NULL, x, n, 0);
    if (err != NULL)
        *err = 15;
}

void _mp_real_cosh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_sc_rs(NULL, res, x, n, 0);
    if (err != NULL)
        *err = 15;
}

void _mp_real_atan_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_at_rs(res, x, n, 1);
    if (err != NULL)
        *err = 15;
}

void _mp_real_atanh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_at_rs(res, x, n, 0);
    if (err != NULL)
        *err = 15;
}
