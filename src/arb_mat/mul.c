/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_mul(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong cutoff;

    /* todo: detect small-integer matrices */
    if (prec <= 2 * FLINT_BITS)
        cutoff = 60;
    else if (prec <= 8 * FLINT_BITS)
        cutoff = 50;
    else
        cutoff = 40;

    if (arb_mat_nrows(A) <= cutoff || arb_mat_ncols(A) <= cutoff ||
        arb_mat_ncols(B) <= cutoff)
    {
        if (flint_get_num_threads() > 1 &&
            ((double) arb_mat_nrows(A) *
             (double) arb_mat_nrows(B) *
             (double) arb_mat_ncols(B) *
             (double) prec > 100000))
        {
            arb_mat_mul_threaded(C, A, B, prec);
        }
        else
        {
            arb_mat_mul_classical(C, A, B, prec);
        }
    }
    else
    {
        arb_mat_mul_block(C, A, B, prec);
    }
}
