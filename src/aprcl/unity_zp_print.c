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
unity_zp_print(const unity_zp f)
{
    flint_printf("p = %wu; exp = %wu\n", f->p, f->exp);
    fmpz_mod_poly_print(f->poly, f->ctx);
    flint_printf("\n");
}

