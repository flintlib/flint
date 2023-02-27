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
unity_zp_coeff_inc(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_ui(f->poly, ind, 1, f->ctx);
        return;
    }

    fmpz_add_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal(f->poly->coeffs + ind, fmpz_mod_ctx_modulus(f->ctx)))
        fmpz_set_ui(f->poly->coeffs + ind, 0);
}

void
unity_zp_coeff_dec(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_si(f->poly, ind, -1, f->ctx);
        return;
    }

    fmpz_sub_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal_si(f->poly->coeffs + ind, -1))
        fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

