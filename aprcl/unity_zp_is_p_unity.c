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

slong
unity_zp_is_p_unity(const unity_zp f)
{
    ulong i, j, m, size, p_pow1, p_pow2;
    slong h;
    fmpz_t n;

    fmpz_init_set(n, f->n);
    fmpz_sub_ui(n, n, 1);
    p_pow1 = n_pow(f->p, f->exp - 1);
    p_pow2 = p_pow1 * f->p;
    h = -1;
    
    size = FLINT_MIN(p_pow2, f->poly->length);
    for (i = 0; i < size; i++)
    {
        if (fmpz_equal_ui(f->poly->coeffs + i, 0) == 0 && h != -1)
        {
            h = -1;
            break;
        }
        if (fmpz_equal_ui(f->poly->coeffs + i, 1))
            h = i;
    }

    /* reduce cyclotomic... */

    fmpz_clear(n);
    return h;
}

