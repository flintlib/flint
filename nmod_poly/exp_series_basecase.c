/*
    Copyright (C) 2011, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_exp_series_basecase(mp_ptr f, mp_srcptr h,
                                    slong hlen, slong n, nmod_t mod)
{
    slong j, k, l;
    mp_ptr a;
    mp_limb_t s;
    int nlimbs;
    TMP_INIT;

    f[0] = 1;

    if (hlen < 2)
    {
        _nmod_vec_zero(f + 1, n - 1);
        return;
    }

    if (n < 2)
        return;

    f[1] = h[1];

    TMP_START;
    a = TMP_ALLOC(FLINT_MIN(n, hlen) * sizeof(mp_limb_t));

    for (k = 1; k < FLINT_MIN(n, hlen); k++)
        a[k] = n_mulmod2_preinv(h[k], k, mod.n, mod.ninv);

    nlimbs = _nmod_vec_dot_bound_limbs(FLINT_MIN(n, hlen), mod);

    for (k = 2; k < n; k++)
    {
        l = FLINT_MIN(k, hlen - 1);
        NMOD_VEC_DOT(s, j, l, a[1 + j], f[k - 1 - j], mod, nlimbs)
        f[k] = n_mulmod2_preinv(s, n_invmod(k, mod.n), mod.n, mod.ninv);
    }

    TMP_END;
}

void
nmod_poly_exp_series_basecase(nmod_poly_t f, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != 0)
    {
        flint_printf("Exception (nmod_poly_exp_series_basecase). Constant term != 0.\n");
        flint_abort();
    }

    if (n <= 1 || hlen <= 1)
    {
        if (n == 0)
            nmod_poly_zero(f);
        else
            nmod_poly_one(f);
        return;
    }

    nmod_poly_fit_length(f, n);
    _nmod_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, f->mod);
    f->length = n;
    _nmod_poly_normalise(f);
}
