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
unity_zpq_coeff_set_fmpz(unity_zpq f, slong i, slong j, const fmpz_t x)
{
    fmpz_mod_poly_set_coeff_fmpz(f->polys[j], i, x, f->ctx);
}

void
unity_zpq_coeff_set_ui(unity_zpq f, slong i, slong j, ulong x)
{
    fmpz_mod_poly_set_coeff_ui(f->polys[j], i, x, f->ctx);
}

