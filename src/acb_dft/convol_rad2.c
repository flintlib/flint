/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

/* assume np >= 2 * n - 1 */
void
acb_dft_convol_pad(acb_ptr fp, acb_ptr gp, acb_srcptr f, acb_srcptr g, slong n, slong np)
{
    slong k;

    if (np < 2 * n - 1)
    {
        flint_throw(FLINT_ERROR, "dft_convol_pad: overlapping padding %wd < 2*%wd-1\n", np, n);
    }

    for (k = 0; k < n; k++)
        acb_set(gp + k, g + k);
    for (; k < np; k++)
        acb_zero(gp + k);

    for (k = 0; k < n; k++)
        acb_set(fp + k, f + k);
    for (k = 1; k < n; k++)
        acb_set(fp + np - k, f + n - k);
    for (k = n; k <= np - n; k++)
        acb_zero(fp + k);

}

void
acb_dft_convol_rad2_precomp(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, const acb_dft_rad2_t rad2, slong prec)
{
    slong np;
    acb_ptr fp, gp;
    np = rad2->n;

    if (len <= 0)
        return;

    fp = _acb_vec_init(np);
    gp = _acb_vec_init(np);

    if (len == np)
    {
        _acb_vec_set(fp, f, len);
        _acb_vec_set(gp, g, len);
    }
    else
        acb_dft_convol_pad(fp, gp, f, g, len, np);

    acb_dft_rad2_precomp_inplace(fp, rad2, prec);
    acb_dft_rad2_precomp_inplace(gp, rad2, prec);

    _acb_vec_kronecker_mul(gp, gp, fp, np, prec);

    acb_dft_inverse_rad2_precomp_inplace(gp, rad2, prec);

    _acb_vec_set(w, gp, len);
    _acb_vec_clear(fp, np);
    _acb_vec_clear(gp, np);
}

void
acb_dft_convol_rad2(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    int e;
    acb_dft_rad2_t dft;
    /* catch power of 2 */
    if (len <= 0)
        return;
    else if ((len & (len - 1)) == 0)
        e = n_clog(len, 2);
    else
        e = n_clog(2 * len - 1, 2);
    acb_dft_rad2_init(dft, e, prec);
    acb_dft_convol_rad2_precomp(w, f, g, len, dft, prec);
    acb_dft_rad2_clear(dft);
}
