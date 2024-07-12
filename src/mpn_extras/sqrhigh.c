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

#include "mpn_extras.h"

static const short flint_mpn_sqrhigh_k_tab[FLINT_MPN_SQRHIGH_K_TAB_SIZE] = {FLINT_MPN_SQRHIGH_K_TAB};

void
_flint_mpn_sqrhigh_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_size_t n)
{
    mp_limb_t cy;
    mp_size_t l;
    slong k;

    if (FLINT_HAVE_SQRHIGH_FUNC(n))
    {
        rp[n - 1] = flint_mpn_sqrhigh_func_tab[n](rp + n, np);
        return;
    }

    if (n < FLINT_MPN_SQRHIGH_K_TAB_SIZE)
        k = flint_mpn_sqrhigh_k_tab[n];
    else
        k = (n + 4) / 2;

    if (k == 0)
    {
        rp[n - 1] = _flint_mpn_sqrhigh_basecase(rp + n, np, n);
        return;
    }

    FLINT_ASSERT(k >= (n + 4) / 2);

    if (k == n)
    {
        flint_mpn_sqr(rp, np, n);
        return;
    }

    l = n - k;

    flint_mpn_sqr(rp + 2 * l, np + l, k);            /* fills rp[2l..2n-1] */
    _flint_mpn_mulhigh_n_mulders_recursive(rp, np, np + k, l);         /* fills rp[l-1..2l-1] */
    /* {rp+n-1,l+1} += 2 * {rp+l-1,l+1} */
    cy = mpn_lshift(rp + l - 1, rp + l - 1, l + 1, 1);
    cy += mpn_add_n(rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
    mpn_add_1(rp + n + l, rp + n + l, k, cy); /* propagate carry */
}

mp_limb_t
_flint_mpn_sqrhigh_mulders(mp_ptr res, mp_srcptr u, mp_size_t n)
{
    mp_ptr tmp, tr, tu;
    mp_limb_t bot;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (3 * (n + 1)));
    tu = tmp;
    tr = tu + (n + 1);
    tu[0] = 0;
    flint_mpn_copyi(tu + 1, u, n);
    _flint_mpn_sqrhigh_mulders_recursive(tr, tu, n + 1);
    flint_mpn_copyi(res, tr + n + 2, n);
    bot = tr[n + 1];

    TMP_END;
    return bot;
}

mp_limb_t
_flint_mpn_sqrhigh_sqr(mp_ptr res, mp_srcptr u, mp_size_t n)
{
    mp_ptr tmp;
    mp_limb_t bot;
    tmp = flint_malloc(sizeof(mp_limb_t) * (2 * n));
    flint_mpn_sqr(tmp, u, n);
    flint_mpn_copyi(res, tmp + n, n);
    bot = tmp[n - 1];
    flint_free(tmp);
    return bot;
}

mp_limb_t
_flint_mpn_sqrhigh(mp_ptr res, mp_srcptr u, mp_size_t n)
{
    if (n <= FLINT_MPN_SQRHIGH_MULDERS_CUTOFF)
        return _flint_mpn_sqrhigh_basecase(res, u, n);
    else if (n <= FLINT_MPN_SQRHIGH_SQR_CUTOFF)
        return _flint_mpn_sqrhigh_mulders(res, u, n);
    else
        return _flint_mpn_sqrhigh_sqr(res, u, n);
}

mp_limb_pair_t _flint_mpn_sqrhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_size_t n)
{
    mp_limb_pair_t ret;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(rp != xp);

    ret.m1 = flint_mpn_sqrhigh(rp, xp, n);

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
