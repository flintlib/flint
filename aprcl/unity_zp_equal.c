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

int
unity_zp_equal(unity_zp f, unity_zp g)
{
    /*
        f and g can be reduced only by modylo x^{p^k} - 1,
        so reduce by cyclotomic polynomial
    */
    _unity_zp_reduce_cyclotomic(f);
    _unity_zp_reduce_cyclotomic(g);

    if (fmpz_mod_poly_equal(f->poly, g->poly) == 0)
        return 0;

    return 1;
}

