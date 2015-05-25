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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void 
fmpz_ec_add(fmpz_t x, fmpz_t y, fmpz_t x1, fmpz_t z1, fmpz_t x2, fmpz_t z2,
            fmpz_t x0, fmpz_t z0)
{
    /* x = x0 * ((x2 - z2) * (x1 + z1) + (x2 + z2) * (x1 - z1))^2 */
    /* z = z0 * ((x2 - z2) * (x1 + z1) + (x2 + z2) * (x1 - z1))^2 */

    fmpz_t u, v, w;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(w);

    fmpz_sub(u, x2, z2);
    fmpz_add(v, x1, z1);
    fmpz_mul(u, u, v);

    fmpz_add(w, x2, z2);
    fmpz_sub(v, x1, z1);
    fmpz_mul(v, w, v);

    fmpz_add(w, u, v);
    fmpz_sub(v, u, v);

    fmpz_mul(w, w, w);
    fmpz_mul(v, v, v);
    
    fmpz_mul(x, z0, w);
    fmpz_mul(z, x0, v);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(w);

}

