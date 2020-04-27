/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_set_fmpz(ca_qqbar_t res, const fmpz_t x)
{
    fmpz_poly_zero(CA_QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(CA_QQBAR_POLY(res), 1, 1);
    fmpz_neg(CA_QQBAR_POLY(res)->coeffs, x);
    acb_set_round_fmpz(CA_QQBAR_ENCLOSURE(res), x, CA_QQBAR_DEFAULT_PREC);
}

