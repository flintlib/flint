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

void
unity_zp_pow_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
{
    slong i;
    unity_zp temp;

    unity_zp_init(temp, f->p, f->exp, f->n);

    unity_zp_set_zero(f);
    unity_zp_coeff_set_ui(f, 0, 1);

    i = fmpz_bits(pow) - 1;

    while (i >= 0)
    {
        unity_zp_sqr(temp, f);
        unity_zp_swap(f, temp);

        if (fmpz_tstbit(pow, i) == 1)
        {
            unity_zp_mul(temp, f, g);
            unity_zp_swap(f, temp);
        }

        i--;
    }

    unity_zp_clear(temp);
}

void
unity_zp_pow_ui(unity_zp f, const unity_zp g, ulong pow)
{
    fmpz_t p;
    fmpz_init_set_ui(p, pow);
    unity_zp_pow_fmpz(f, g, p);
    fmpz_clear(p);
}

