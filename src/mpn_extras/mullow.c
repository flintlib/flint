/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpn_extras.h"

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
mp_limb_t flint_mpn_mullow_1(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_2(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_3(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_4(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_5(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_6(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_7(mp_ptr, mp_srcptr, mp_srcptr);
mp_limb_t flint_mpn_mullow_8(mp_ptr, mp_srcptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mullow_func_tab[] =
{
    NULL,
    flint_mpn_mullow_1,
    flint_mpn_mullow_2,
    flint_mpn_mullow_3,
    flint_mpn_mullow_4,
    flint_mpn_mullow_5,
    flint_mpn_mullow_6,
    flint_mpn_mullow_7,
    flint_mpn_mullow_8
};
#else
const flint_mpn_mul_func_t flint_mpn_mullow_func_tab[] = { NULL };
#endif

#if !FLINT_HAVE_NATIVE_mpn_mullow_basecase
mp_limb_t
flint_mpn_mullow_basecase(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_limb_t ret;
    mp_size_t ix;

    ret = mpn_mul_1(rp, xp, n, yp[0]);

    for (ix = 1; ix < n; ix++)
    {
        ret += mpn_addmul_1(rp + ix, xp, n - ix, yp[ix]);
        ret += xp[n - ix] * yp[ix];
    }

    return ret;
}
#endif

void
_flint_mpn_mullow_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
    mp_size_t l;
    slong k;

    if (FLINT_HAVE_MULLOW_FUNC(n))
    {
        flint_mpn_mullow_func_tab[n](rp, np, mp);
        return;
    }

    if (n < FLINT_MPN_MULHIGH_K_TAB_SIZE)
        k = flint_mpn_mulhigh_k_tab[n];
    else
        k = 3 * (n / 4);

    if (k == 0)
    {
        flint_mpn_mullow_basecase(rp, np, mp, n);
        return;
    }

    if (k == n)
    {
        flint_mpn_mul_n(rp, np, mp, n);
        return;
    }

    FLINT_ASSERT(k >= (n + 1) / 2);

    l = n - k;

    flint_mpn_mul_n(rp, np, mp, k);
    _flint_mpn_mullow_n_mulders_recursive(rp + n, np, mp + k, l);
    mpn_add_n(rp + k, rp + k, rp + n, l);
    _flint_mpn_mullow_n_mulders_recursive(rp + n, np + k, mp, l);
    mpn_add_n(rp + k, rp + k, rp + n, l);
}

mp_limb_t
_flint_mpn_mullow_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_ptr tmp, tr, tu, tv;
    mp_limb_t bot;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (4 * (n + 1)));
    tu = tmp;
    tv = tu + (n + 1);
    tr = tv + (n + 1);
    tu[n] = 0;
    tv[n] = 0;
    flint_mpn_copyi(tu, u, n);
    flint_mpn_copyi(tv, v, n);
    _flint_mpn_mullow_n_mulders_recursive(tr, tu, tv, n + 1);
    flint_mpn_copyi(res, tr, n);
    bot = tr[n];

    TMP_END;
    return bot;
}

mp_limb_t
_flint_mpn_mullow_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_ptr tmp;
    mp_limb_t cy;
    tmp = flint_malloc(sizeof(mp_limb_t) * (2 * n));
    flint_mpn_mul_n(tmp, u, v, n);
    memcpy(res, tmp, sizeof(mp_limb_t) * n);
    cy = tmp[n];
    flint_free(tmp);
    return cy;
}

mp_limb_t
_flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    if (n <= FLINT_MPN_MULHIGH_MULDERS_CUTOFF)
        return flint_mpn_mullow_basecase(res, u, v, n);
    else if (n <= FLINT_MPN_MULHIGH_MUL_CUTOFF)
        return _flint_mpn_mullow_n_mulders(res, u, v, n);
    else
        return _flint_mpn_mullow_n_mul(res, u, v, n);
}
