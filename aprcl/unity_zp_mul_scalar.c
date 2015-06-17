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
unity_zp_mul_scalar_fmpz(unity_zp f, const unity_zp g, const fmpz_t s)
{
    fmpz_mod_poly_scalar_mul_fmpz(f->poly, g->poly, s);
}

void
unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)
{
    fmpz_mod_poly_scalar_mul_ui(f->poly, g->poly, s);
}
