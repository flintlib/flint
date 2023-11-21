/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_phi(qqbar_t res)
{
    fmpz_poly_zero(QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 2, 1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 1, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 0, -1);

    arb_sqrt_ui(acb_realref(QQBAR_ENCLOSURE(res)), 5, QQBAR_DEFAULT_PREC);
    arb_add_ui(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), 1, QQBAR_DEFAULT_PREC);
    arb_mul_2exp_si(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), -1);
    arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
}

