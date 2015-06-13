/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"

void _unity_zp_reduce_cyclotomic(unity_zp f)
{
    ulong i, j, ppow, cycl_pow;
    fmpz_t coeff;

    ppow = n_pow(f->p, f->exp - 1);
    cycl_pow = (f->p - 1) * ppow;
    fmpz_init(coeff);

    for (i = f->poly->length - 1; i >= cycl_pow; i--)
    {
        fmpz_set(coeff, f->poly->coeffs + i);
        if (fmpz_is_zero(coeff))
            continue;

        for (j = 0; j < f->p - 1; j++)
        {
            ulong ind = j * ppow;
            fmpz_sub(f->poly->coeffs + ind
                , f->poly->coeffs + ind, coeff);
            if (fmpz_cmp(f->poly->coeffs + ind, 0) < 0)
                fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind, f->n);
        }
    } 

    fmpz_clear(coeff);
}

void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)
{
    unity_zp_copy(f, g);
    _unity_zp_reduce_cyclotomic(f);
}

