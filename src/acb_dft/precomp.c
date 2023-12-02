/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

void
_acb_dft_precomp_init(acb_dft_pre_t pre, slong dv, acb_ptr z, slong dz, slong len, slong prec)
{
    pre->n = len;
    if (len <= 1)
    {
        pre->type = DFT_NAIVE;
        _acb_dft_naive_init(pre->t.naive, dv, z, dz, len, prec);
    }
    else if (n_is_prime(len))
    {
        if (len < 100)
        {
            pre->type = DFT_NAIVE;
            _acb_dft_naive_init(pre->t.naive, dv, z, dz, len, prec);
        }
        else
        {
            pre->type = DFT_CONV;
            /* FIXME: do not recompute the same bluestein
             * scheme if needed several times */
            _acb_dft_bluestein_init(pre->t.bluestein, dv, len, prec);
        }
    }
    else
    {
        n_factor_t fac;

        n_factor_init(&fac);
        n_factor(&fac, len, 1);

        if (fac.num == 1)
        {
            /* TODO: could be p^e, or 2^e, but with dv shift */
            if (fac.p[0] == 2)
            {
                pre->type = DFT_RAD2;
                _acb_dft_rad2_init(pre->t.rad2, dv, fac.exp[0], prec);
            }
            else
            {
                pre->type = DFT_CYC;
                _acb_dft_cyc_init_z_fac(pre->t.cyc, fac, dv, z, dz, len, prec);
            }
        }
        else
        {
            pre->type = DFT_CRT;
            _acb_dft_crt_init(pre->t.crt, dv, len, prec);
        }
    }
}

void
acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)
{
    _acb_dft_precomp_init(pre, 1, NULL, 0, len, prec);
}

void
acb_dft_precomp_clear(acb_dft_pre_t pre)
{
    switch (pre->type)
    {
        case DFT_NAIVE:
            acb_dft_naive_clear(pre->t.naive);
            break;
        case DFT_CYC:
            acb_dft_cyc_clear(pre->t.cyc);
            break;
        case DFT_PROD:
            acb_dft_prod_clear(pre->t.prod);
            break;
        case DFT_CRT:
            acb_dft_crt_clear(pre->t.crt);
            break;
        case DFT_RAD2:
            acb_dft_rad2_clear(pre->t.rad2);
            break;
        case DFT_CONV:
            acb_dft_bluestein_clear(pre->t.bluestein);
            break;
        default:
            flint_throw(FLINT_ERROR, "acb_dft_clear: unknown strategy code %i\n", pre->type);
    }
}

void
acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
{
    switch (pre->type)
    {
        case DFT_NAIVE:
            acb_dft_naive_precomp(w, v, pre->t.naive, prec);
            break;
        case DFT_CYC:
            acb_dft_cyc_precomp(w, v, pre->t.cyc, prec);
            break;
        case DFT_PROD:
            acb_dft_prod_precomp(w, v, pre->t.prod, prec);
            break;
        case DFT_CRT:
            acb_dft_crt_precomp(w, v, pre->t.crt, prec);
            break;
        case DFT_RAD2:
            acb_dft_rad2_precomp(w, v, pre->t.rad2, prec);
            break;
        case DFT_CONV:
            acb_dft_bluestein_precomp(w, v, pre->t.bluestein, prec);
            break;
        default:
            flint_throw(FLINT_ERROR, "acb_dft_precomp: unknown strategy code %i\n", pre->type);
    }
}

void
acb_dft(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dft_pre_t t;
    acb_dft_precomp_init(t, len, prec);
    acb_dft_precomp(w, v, t, prec);
    acb_dft_precomp_clear(t);
}

void
acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
{
    slong k;
    acb_dft_precomp(w, v, pre, prec);
    for (k = 1; 2 * k < pre->n; k++)
        acb_swap(w + k, w + pre->n - k);
    _acb_vec_scalar_div_ui(w, w, pre->n, pre->n, prec);
}
void
acb_dft_inverse(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    slong k;
    acb_dft(w, v, len, prec);
    for (k = 1; 2 * k < len; k++)
        acb_swap(w + k, w + len - k);
    _acb_vec_scalar_div_ui(w, w, len, len, prec);
}
