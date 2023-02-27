/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

void
fmpr_printd(const fmpr_t x, slong digits)
{
    mpfr_t t;
    mpfr_init2(t, digits * 3.33 + 10);
    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);
    fmpr_get_mpfr(t, x, MPFR_RNDN);
    mpfr_printf("%.*Rg", FLINT_MAX(digits, 1), t);
    mpfr_clear(t);
}
