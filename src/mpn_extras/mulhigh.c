/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    Uses code adopted from the GNU MPFR Library.

        Copyright 2005-2024 Free Software Foundation, Inc.
        Contributed by the AriC and Caramba projects, INRIA.

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpn_extras.h"

#if TUNE_PROGRAM
short flint_mpn_mulhigh_k_tab[FLINT_MPN_MULHIGH_K_TAB_SIZE];
#else
const short flint_mpn_mulhigh_k_tab[FLINT_MPN_MULHIGH_K_TAB_SIZE] = {FLINT_MPN_MULHIGH_K_TAB};
#endif

void
_flint_mpn_mulhigh_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
    mp_limb_t cy;
    mp_size_t l;
    slong k;

    if (FLINT_HAVE_MULHIGH_FUNC(n))
    {
        rp[n - 1] = flint_mpn_mulhigh_func_tab[n](rp + n, np, mp);
        return;
    }

    if (n < FLINT_MPN_MULHIGH_K_TAB_SIZE)
        k = flint_mpn_mulhigh_k_tab[n];
    else
        k = 3 * (n / 4);

    if (k == 0)
    {
        rp[n - 1] = _flint_mpn_mulhigh_basecase(rp + n, np, mp, n);
        return;
    }

    FLINT_ASSERT(k >= (n + 4) / 2);

    if (k == n)
    {
        flint_mpn_mul_n(rp, np, mp, n);
        return;
    }

    l = n - k;
    flint_mpn_mul_n(rp + 2 * l, np + l, mp + l, k);
    _flint_mpn_mulhigh_n_mulders_recursive(rp, np + k, mp, l);
    cy = mpn_add_n(rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
    _flint_mpn_mulhigh_n_mulders_recursive(rp, np, mp + k, l);
    cy += mpn_add_n(rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
    MPN_INCR_U(rp + n + l, k, cy);
}

mp_limb_t
_flint_mpn_mulhigh_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_ptr tmp, tr, tu, tv;
    mp_limb_t bot;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (4 * (n + 1)));
    tu = tmp;
    tv = tu + (n + 1);
    tr = tv + (n + 1);
    tu[0] = 0;
    tv[0] = 0;
    flint_mpn_copyi(tu + 1, u, n);
    flint_mpn_copyi(tv + 1, v, n);
    _flint_mpn_mulhigh_n_mulders_recursive(tr, tu, tv, n + 1);
    flint_mpn_copyi(res, tr + n + 2, n);
    bot = tr[n + 1];

    TMP_END;
    return bot;
}

mp_limb_t
_flint_mpn_mulhigh_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_ptr tmp;
    mp_limb_t bot;
    tmp = flint_malloc(sizeof(mp_limb_t) * (2 * n));
    flint_mpn_mul_n(tmp, u, v, n);
    memcpy(res, tmp + n, sizeof(mp_limb_t) * n);
    bot = tmp[n - 1];
    flint_free(tmp);
    return bot;
}

#if !TUNE_PROGRAM
mp_limb_t
_flint_mpn_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    if (n <= FLINT_MPN_MULHIGH_MULDERS_CUTOFF)
        return _flint_mpn_mulhigh_basecase(res, u, v, n);
    else if (n <= FLINT_MPN_MULHIGH_MUL_CUTOFF)
        return _flint_mpn_mulhigh_n_mulders(res, u, v, n);
    else
        return _flint_mpn_mulhigh_n_mul(res, u, v, n);
}

mp_limb_pair_t _flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    mp_limb_pair_t ret;

    FLINT_ASSERT(rp != xp && rp != yp);

    ret.m1 = flint_mpn_mulhigh_n(rp, xp, yp, n);

    if (rp[n - 1] >> (FLINT_BITS - 1))
    {
        ret.m2 = 0;
    }
    else
    {
        ret.m2 = 1;
        mpn_lshift(rp, rp, n, 1);
        rp[0] |= (ret.m1 >> (FLINT_BITS - 1));
        ret.m1 <<= 1;
    }

    return ret;
}
#endif
