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
ca_qqbar_set_fmpq(ca_qqbar_t res, const fmpq_t x)
{
    fmpz_poly_zero(CA_QQBAR_POLY(res));
    fmpz_poly_set_coeff_fmpz(CA_QQBAR_POLY(res), 1, fmpq_denref(x));
    fmpz_neg(CA_QQBAR_COEFFS(res), fmpq_numref(x));
    acb_set_fmpq(CA_QQBAR_ENCLOSURE(res), x, CA_QQBAR_DEFAULT_PREC);
}

