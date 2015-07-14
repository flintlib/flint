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

/* Init Montgomery form structure */
void
unity_zp_mont_init(unity_zp_mont f, ulong p, ulong exp, const fmpz_t n, const fmpz_t ninv)
{
    f->p = p;
    f->exp = exp;
    fmpz_init_set(f->n, n);
    fmpz_init(f->nr);

    fmpz_init_set(f->ninv, ninv);
    f->r = fmpz_bits(n);

    fmpz_poly_init(f->poly);
}

void
unity_zp_mont_clear(unity_zp_mont f)
{
    fmpz_clear(f->nr);
    fmpz_clear(f->n);
    fmpz_clear(f->ninv);
    fmpz_poly_clear(f->poly);
}

