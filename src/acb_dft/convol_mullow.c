/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"
#include "acb_poly.h"

void
acb_dft_convol_mullow(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    /* TODO: should probably use (acb_struct *) arrays */
    acb_ptr gg, ww;
    if (len == 0)
        return;
    gg = _acb_vec_init(2 * len - 1);
    ww = _acb_vec_init(2 * len - 1);
    _acb_vec_set(gg, g, len);
    _acb_vec_set(gg + len, g, len - 1);
    _acb_poly_mullow(ww, f, len, gg, 2 * len - 1, 2 * len - 1, prec);
    _acb_vec_set(w, ww + len, len - 1);
    acb_set(w + len - 1, ww + len - 1);
    _acb_vec_clear(gg, 2 * len - 1);
    _acb_vec_clear(ww, 2 * len - 1);
}
