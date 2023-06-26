/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_spd_eig_lbound_arf(arf_t b, const arb_mat_t mat, slong prec)
{
    arb_poly_t poly;
    arb_t x;

    arb_poly_init(poly);
    arb_init(x);

    arb_mat_charpoly(poly, mat, prec);
    arb_div(x, arb_poly_get_coeff_ptr(poly, 0),
            arb_poly_get_coeff_ptr(poly, 1), prec);
    arb_neg(x, x);
    arb_get_lbound_arf(b, x, prec);

    arb_poly_clear(poly);
    arb_clear(x);
}
