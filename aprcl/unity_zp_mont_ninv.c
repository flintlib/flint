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
unity_zp_mont_ninv(fmpz_t ninv, const unity_zp_mont f)
{
    fmpz_t d, n, r, rinv;

    fmpz_init(d);
    fmpz_init(n);
    fmpz_init(r);
    fmpz_init(rinv);

    fmpz_set(n, f->n);
    fmpz_set_ui(d, 1);
    fmpz_setbit(r, fmpz_bits(n));
    fmpz_neg(n, n);


    fmpz_xgcd(d, rinv, ninv, r, n);
    fmpz_mod(ninv, ninv, r);

    fmpz_clear(d);
    fmpz_clear(n);
    fmpz_clear(r);
    fmpz_clear(rinv);
}

