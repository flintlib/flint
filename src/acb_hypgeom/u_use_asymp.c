/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

int
acb_hypgeom_u_use_asymp(const acb_t z, slong prec)
{
    double x, y;

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0))
    {
        return 0;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        return 1;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    return sqrt(x * x + y * y) > prec * 0.69314718055994530942;
}

