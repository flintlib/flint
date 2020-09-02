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
_unity_zp_reduce_cyclotomic_divmod(unity_zp f)
{
    ulong i, j, ppow1, ppow2, cycl_pow;

    ppow2 = n_pow(f->p, f->exp - 1);
    ppow1 = ppow2 * f->p;
    cycl_pow = (f->p - 1) * ppow2;

    for (i = f->poly->length - 1; i >= ppow1; i--)
    {
        fmpz_add(f->poly->coeffs + i - ppow1,
                f->poly->coeffs + i - ppow1, f->poly->coeffs + i);

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    for (i = f->poly->length - 1; i >= cycl_pow; i--)
    {
        if (fmpz_is_zero(f->poly->coeffs + i))
            continue;

        for (j = 0; j < f->p - 1; j++)
        {
            ulong ind = i - cycl_pow + j * ppow2;
            fmpz_sub(f->poly->coeffs + ind,
                    f->poly->coeffs + ind, f->poly->coeffs + i);
        }

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _fmpz_mod_poly_normalise(f->poly);
    _fmpz_vec_scalar_mod_fmpz(f->poly->coeffs,
               f->poly->coeffs, f->poly->length, fmpz_mod_ctx_modulus(f->ctx));
    _fmpz_mod_poly_normalise(f->poly);
}

void
_unity_zp_reduce_cyclotomic(unity_zp f)
{
    ulong i, j, ppow, cycl_pow;

    if (f->poly->length == 0)
        return;

    ppow = n_pow(f->p, f->exp - 1);
    cycl_pow = (f->p - 1) * ppow;

    for (i = f->poly->length - 1; i >= cycl_pow; i--)
    {
        if (fmpz_is_zero(f->poly->coeffs + i))
            continue;

        for (j = 0; j < f->p - 1; j++)
        {
            ulong ind = i - cycl_pow + j * ppow;
            fmpz_sub(f->poly->coeffs + ind,
                    f->poly->coeffs + ind, f->poly->coeffs + i);

            if (fmpz_cmp_ui(f->poly->coeffs + ind, 0) < 0)
                fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
        }

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _fmpz_mod_poly_normalise(f->poly);
}

void
unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)
{
    unity_zp_copy(f, g);
    _unity_zp_reduce_cyclotomic(f);
}
