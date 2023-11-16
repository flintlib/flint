/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "acb_dft.h"
#include "acb_modular.h"

/* z[k] = z^(k^2), z a 2n-th root of unity */
static void
_acb_vec_bluestein_factors(acb_ptr z, slong n, slong prec)
{
    /* this function is used mostly with prime-power n
     * so the set of squares has index 2 only.
     * computing an addition sequence does not considerably improve things */
    if (n < 30)
    {
        slong k, k2;
        acb_ptr z2n;
        nmod_t n2;

        z2n = _acb_vec_init(2 * n);
        _acb_vec_unit_roots(z2n, -2 * n, 2 * n, prec);
        nmod_init(&n2, FLINT_MAX(2 * n, 1));

        for (k = 0, k2 = 0; k < n; k++)
        {
            acb_set(z + k, z2n + k2);
            k2 = nmod_add(k2, 2 * k + 1, n2);
        }

        _acb_vec_clear(z2n, 2 * n);
    }
    else
    {
        nmod_t n2;
        slong k, k2, dk;
        slong * v, * s;
        acb_ptr t;
        s = flint_malloc(n * sizeof(slong));
        v = flint_malloc((n + 1)* sizeof(slong));
        t = _acb_vec_init(n + 1);
        nmod_init(&n2, 2 * n);

        for (k = 0; k < n; k++)
            v[k] = 0;
        for (k = 0, k2 = 0, dk = 1; k < n; k++)
        {
            s[k] = k2;
            if (k2 < n)
                v[k2] = -1;
            else
                v[2 * n - k2] = -1;

            k2 = nmod_add(k2, dk, n2);
            dk = nmod_add(dk, 2, n2);
        }
        acb_modular_fill_addseq(v, n);

        acb_one(t + 0);
        acb_unit_root(t + 1, 2 * n, prec);
        acb_conj(t + 1, t + 1);
        acb_set_si(t + n, -1);
        for (k = 2; k < n; k++)
            if (v[k])
                acb_mul(t + k, t + v[k], t + k - v[k], prec);
        for (k = 0; k < n; k++)
        {
            if (s[k] <= n)
                acb_set(z + k, t + s[k]);
            else
                acb_conj(z + k, t + 2 * n - s[k]);
        }
        _acb_vec_clear(t, n + 1);
        flint_free(s);
        flint_free(v);
    }
}

void
_acb_dft_bluestein_init(acb_dft_bluestein_t t, slong dv, slong n, slong prec)
{
    acb_ptr z, g;
    slong k, n2;
    int e;

    t->n = n;
    t->dv = dv;

    if (n == 0)
        return;

    e = n_clog(2 * n - 1, 2);

    if (DFT_VERB)
        flint_printf("dft_bluestein: init z[2^%i]\n", e);

    acb_dft_rad2_init(t->rad2, e, prec);

    t->z = z = _acb_vec_init(n);
    _acb_vec_bluestein_factors(t->z, n, prec);

    n2 = t->rad2->n;
    t->g = g = _acb_vec_init(n2);
    acb_one(g + 0);
    for (k = 1; k < n; k++)
    {
        acb_conj(g + k, z + k);
        acb_set(g + n2 - k, g + k);
    }
    acb_dft_rad2_precomp_inplace(g, t->rad2, prec);
}

void
acb_dft_bluestein_precomp(acb_ptr w, acb_srcptr v, const acb_dft_bluestein_t t, slong prec)
{
    slong n = t->n, np = t->rad2->n, dv = t->dv;
    acb_ptr fp;

    if (n == 0)
        return;

    fp = _acb_vec_init(np);
    _acb_vec_kronecker_mul_step(fp, t->z, v, dv, n, prec);

    acb_dft_rad2_precomp_inplace(fp, t->rad2, prec);
    _acb_vec_kronecker_mul(fp, t->g, fp, np, prec);

    acb_dft_inverse_rad2_precomp_inplace(fp, t->rad2, prec);

    _acb_vec_kronecker_mul(w, t->z, fp, n, prec);

    _acb_vec_clear(fp, n);
}

void
acb_dft_bluestein(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dft_bluestein_t t;
    acb_dft_bluestein_init(t, len, prec);
    acb_dft_bluestein_precomp(w, v, t, prec);
    acb_dft_bluestein_clear(t);
}
