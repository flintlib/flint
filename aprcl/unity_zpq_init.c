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

void unity_zpq_init(unity_zpq value, ulong q, ulong p, const fmpz_t n)
{
    int i;

    value->p = p;
    value->q = q;
    fmpz_init_set(value->n, n);
    value->polys = (fmpz_mod_poly_t *) flint_malloc(p * sizeof(fmpz_mod_poly_t));
    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_init(value->polys[i], n);
    }
}

void unity_zpq_clear(unity_zpq value)
{
    int i;
    
    for (i = 0; i < value->p; i++)
    {
        fmpz_mod_poly_clear(value->polys[i]);
    }
    value->p = 0;
    value->q = 0;
    fmpz_clear(value->n);
    flint_free(value->polys);
}

