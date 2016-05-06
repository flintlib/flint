/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zp_mul_scalar_fmpz(unity_zp f, const unity_zp g, const fmpz_t s)
{
    fmpz_mod_poly_scalar_mul_fmpz(f->poly, g->poly, s);
}

void
unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)
{
    fmpz_mod_poly_scalar_mul_ui(f->poly, g->poly, s);
}
