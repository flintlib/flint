/*
    Copyright (C) 2017 Fredrik Johansson
    Copyright (C) 2017 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_csc_pi(arb_t res, const arb_t x, slong prec)
{
    arb_sin_pi(res, x, prec + 4);
    arb_inv(res, res, prec);
}

