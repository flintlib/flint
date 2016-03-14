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
unity_zpq_coeff_set_fmpz(unity_zpq f, ulong i, ulong j, const fmpz_t x)
{
    fmpz_mod_poly_set_coeff_fmpz(f->polys[j], i, x);
}

void
unity_zpq_coeff_set_ui(unity_zpq f, ulong i, ulong j, ulong x)
{
    fmpz_mod_poly_set_coeff_ui(f->polys[j], i, x);
}

