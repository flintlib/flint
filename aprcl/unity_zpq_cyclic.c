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

ulong unity_zpq_punity(const unity_zpq f)
{
    ulong i, is_punity;
    fmpz_t coeff;
    fmpz_init(coeff);
    is_punity = 0;
    for (i = 1; i < f->p; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(coeff, f, i);
        if (fmpz_equal_ui(coeff, 1) == 1)
        {
            if (is_unity != 0)
                return 0;
            is_unity = i;
        }
    }
    fmpz_clear(coeff);
}

int unity_zpq_is_punity(const unity_zpq f)
{
    if (unity_zpq_puinty(f) != 0)
        return 1;
    return 0;
}

int unity_zpq_is_punity_generator(const unity_zpq f)
{
    ulong upow = unity_zpq_punity(f)
    if (upow != 0)
        if (n_gcd(upow, f->p) == 1)
            return 1;
    return 0;
}

