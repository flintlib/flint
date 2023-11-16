/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_sum_fmpq_arb(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
{
    if (N <= 2 || (prec <= 1024 && N <= 8) || (prec <= 4096 && N <= 4))
        arb_hypgeom_sum_fmpq_arb_forward(res, a, alen, b, blen, z, reciprocal, N, prec);
    else if (prec >= 8192 && arb_bits(z) <= 0.001 * prec)
        arb_hypgeom_sum_fmpq_arb_bs(res, a, alen, b, blen, z, reciprocal, N, prec);
    else
        arb_hypgeom_sum_fmpq_arb_rs(res, a, alen, b, blen, z, reciprocal, N, prec);
}
