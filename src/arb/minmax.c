/*
    Copyright (C) 2023 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_minmax(arb_t z1, arb_t z2, const arb_t x, const arb_t y, slong prec)
{
    arf_t xleft, xright, yleft, yright, xr, yr;

    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
    {
        arb_indeterminate(z1);
	arb_indeterminate(z2);
        return;
    }

    if (!arb_is_finite(x) || !arb_is_finite(y))
    {
	if (z1 != x && z1 != y)
	{
	    arb_min(z1, x, y, prec);
	    arb_max(z2, x, y, prec);
	}
	else
	{
	    arb_t t;
	    arb_init(t);
	    arb_min(t, x, y, prec);
	    arb_max(z2, x, y, prec);
	    arb_swap(z1, t);
	    arb_clear(t);
	}
	return;
    }

    arf_init(xleft);
    arf_init(xright);
    arf_init(yleft);
    arf_init(yright);

    arf_init_set_mag_shallow(xr, arb_radref(x));
    arf_init_set_mag_shallow(yr, arb_radref(y));

    arf_sub(xleft, arb_midref(x), xr, prec, ARF_RND_FLOOR);
    arf_sub(yleft, arb_midref(y), yr, prec, ARF_RND_FLOOR);
    arf_add(xright, arb_midref(x), xr, prec, ARF_RND_CEIL);
    arf_add(yright, arb_midref(y), yr, prec, ARF_RND_CEIL);

    if (arf_cmp(xleft, yleft) < 0)
    {
	/* xleft < yleft */
	if (arf_cmp(xright, yright) < 0)
	{
	    /* xright < yright */
	    arb_set_interval_arf(z1, xleft, xright, prec);
	    arb_set_interval_arf(z2, yleft, yright, prec);
	} else {
	    /* xright >= yright */
	    arb_set_interval_arf(z1, xleft, yright, prec);
	    arb_set_interval_arf(z2, yleft, xright, prec);
	}
    } else {
	/* xleft >= yleft */
	if (arf_cmp(xright, yright) < 0)
	{
	    /* xright < yright */
	    arb_set_interval_arf(z1, yleft, xright, prec);
	    arb_set_interval_arf(z2, xleft, yright, prec);
	} else {
	    /* xright >= yright */
	    arb_set_interval_arf(z1, yleft, yright, prec);
	    arb_set_interval_arf(z2, xleft, xright, prec);
	}
    }

    arf_clear(xleft);
    arf_clear(xright);
    arf_clear(yleft);
    arf_clear(yright);
}
