/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
arb_mat_pos_lambda(arb_t lambda, const arb_mat_t mat, slong prec)
{
    arb_poly_t poly;

    arb_poly_init(poly);

    arb_mat_charpoly(poly, mat, prec);
    arb_div(lambda, arb_poly_get_coeff_ptr(poly, 0),
            arb_poly_get_coeff_ptr(poly, 1), prec);
    arb_neg(lambda, lambda);

    arb_poly_clear(poly);
}
