/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x)
{
    fmpz_mod_poly_set_coeff_fmpz(f->poly, ind, x, f->ctx);
}

void
unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x)
{
    fmpz_mod_poly_set_coeff_ui(f->poly, ind, x, f->ctx);
}

