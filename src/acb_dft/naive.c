/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

/* all roots are already computed, w != v */
void
_acb_dft_naive(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong len, slong prec)
{
    slong i, j;
    acb_ptr wi, v1 = NULL;
    acb_srcptr vj;

    if (w == v)
    {
        flint_throw(FLINT_ERROR, "\n_acb_dft_naive: does not accept aliasing\n");
    }

    for (i = 0, wi = w; i < len; i++, wi++)
    {
        acb_zero(wi);
        for (j = 0, vj = v; j < len; j++, vj += dv)
            acb_addmul(wi, vj, z + dz * (i * j % len), prec);
    }

    if (v1)
        _acb_vec_clear(v1, len);
}

void
acb_dft_naive_precomp(acb_ptr w, acb_srcptr v, const acb_dft_naive_t pol, slong prec)
{
    if (v == w)
    {
        acb_ptr v1 = _acb_vec_init(pol->n);
        _acb_vec_set(v1, v, pol->n);
        _acb_dft_naive(w, v1, pol->dv, pol->z, pol->dz, pol->n, prec);
        _acb_vec_clear(v1, pol->n);
    }
    else
        _acb_dft_naive(w, v, pol->dv, pol->z, pol->dz, pol->n, prec);
}

void
acb_dft_naive(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_ptr z, v1 = NULL;

    z = _acb_vec_init(len);
    _acb_vec_unit_roots(z, -len, len, prec);
    if (w == v)
    {
        v1 = _acb_vec_init(len);
        _acb_vec_set(v1, v, len);
        v =  v1;
    }

    _acb_dft_naive(w, v, 1, z, 1, len, prec);

    if (v1)
        _acb_vec_clear(v1, len);

    _acb_vec_clear(z, len);
}

void
_acb_dft_naive_init(acb_dft_naive_t pol, slong dv, acb_ptr z, slong dz, slong len, slong prec)
{
    pol->n = len;
    pol->dv = dv;

    if (z == NULL)
    {
        if (DFT_VERB)
            flint_printf("dft_naive: init z[%wd]\n",len);
        pol->z = _acb_vec_init(len);
        _acb_vec_unit_roots(pol->z, -len, len, prec);
        pol->dz = 1;
        pol->zclear = 1;
    }
    else
    {
        pol->z = z;
        pol->dz = dz;
        pol->zclear = 0;
    }
}
