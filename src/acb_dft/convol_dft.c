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
acb_dft_convol_dft_precomp(acb_ptr w, acb_srcptr f, acb_srcptr g, const acb_dft_pre_t pre, slong prec)
{
    acb_ptr fp, gp;

    fp = _acb_vec_init(pre->n);
    gp = _acb_vec_init(pre->n);

    acb_dft_precomp(fp, f, pre, prec);
    acb_dft_precomp(gp, g, pre, prec);

    _acb_vec_kronecker_mul(gp, gp, fp, pre->n, prec);

    acb_dft_inverse_precomp(w, gp, pre, prec);

    _acb_vec_clear(fp, pre->n);
    _acb_vec_clear(gp, pre->n);
}

void
acb_dft_convol_dft(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    acb_dft_pre_t pre;
    acb_dft_precomp_init(pre, len, prec);
    acb_dft_convol_dft_precomp(w, f, g, pre, prec);
    acb_dft_precomp_clear(pre);
}
