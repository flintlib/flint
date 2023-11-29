/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

/* swap each element with one with bit-reversed index */
void
acb_dft_rad2_reorder(acb_ptr v, slong n)
{
    slong i, j, k, l;

    for (i = 0, l = n>>1; i < l; i++)
    {
        /* j = bit reversal of i */
        for (k = 1, j = 0; k < n; k <<= 1)
        {
            j <<= 1;
            if (i & k)
                j |= 1;
        }
        if (i < j)
            acb_swap(v + i, v + j);
        else if (i > j)
            acb_swap(v + n - 1 - i, v + n - 1 - j);
        i++, j |= l;
        acb_swap(v + i, v + j);
    }

}

void
_acb_dft_rad2_init(acb_dft_rad2_t t, slong dv, int e, slong prec)
{
    if (e < 0)
    {
        flint_throw(FLINT_ERROR, "acb_dft_rad2_init: need e >= 0");
    }
    t->e = e;
    t->n = WORD(1) << e;
    t->dv = dv;
    t->nz = t->n >> 1;
    t->z = _acb_vec_init(t->nz);
    /* set n/2 roots of order n */
    _acb_vec_unit_roots(t->z, -t->n, t->nz, prec);
}

/* remark: can use same rad2 with smaller power of 2 */
void
acb_dft_rad2_precomp_inplace(acb_ptr v, const acb_dft_rad2_t rad2, slong prec)
{
    if (flint_get_num_threads() > 1 && rad2-> e > 9)
    {
        acb_dft_rad2_precomp_inplace_threaded(v, rad2, prec);
    }
    else
    {
        slong j, k, l;
        slong n = rad2->n, nz = rad2->nz;
        acb_ptr p, vend = v + n, w = rad2->z;
        acb_t tmp;
        acb_init(tmp);

        acb_dft_rad2_reorder(v, n);

        for (k = 1, l = nz; k < n; k <<= 1, l >>= 1)
            for (p = v; p < vend; p += k)
                for (j = 0; j < nz; j += l, p++)
                {
                    acb_mul(tmp, p + k, w + j, prec);
                    acb_sub(p + k, p + 0, tmp, prec);
                    acb_add(p + 0, p + 0, tmp, prec);
                }

        acb_clear(tmp);
    }
}

void
acb_dft_inverse_rad2_precomp_inplace(acb_ptr v, const acb_dft_rad2_t rad2, slong prec)
{
    slong k, n = rad2->n;
    acb_dft_rad2_precomp_inplace(v, rad2, prec);
    _acb_vec_scalar_mul_2exp_si(v, v, n, - rad2->e);
    for (k = 1; k < n / 2; k++)
        acb_swap(v + k, v + n - k);
}

void
acb_dft_rad2_inplace(acb_ptr v, int e, slong prec)
{
    acb_dft_rad2_t rad2;
    acb_dft_rad2_init(rad2, e, prec);
    acb_dft_rad2_precomp_inplace(v, rad2, prec);
    acb_dft_rad2_clear(rad2);
}

void
acb_dft_rad2_precomp(acb_ptr w, acb_srcptr v, const acb_dft_rad2_t rad2, slong prec)
{
    slong k;
    for (k = 0; k < rad2->n; k++, v += rad2->dv)
        acb_set(w + k, v + 0);
    acb_dft_rad2_precomp_inplace(w, rad2, prec);
}

void
acb_dft_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)
{
    acb_dft_rad2_t rad2;
    acb_dft_rad2_init(rad2, e, prec);
    acb_dft_rad2_precomp(w, v, rad2, prec);
    acb_dft_rad2_clear(rad2);
}
