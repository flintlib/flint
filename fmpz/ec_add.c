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

/* P (x : z) = P1 (x1 : z1) + P2 (x2 : z2) where P0 (x0 : zo) is P - Q */

/* 
    Coordinates of P : 

        x = 4 * z0 * (x1 * x2 - z1 * z2)^2 
        z = 4 * x0 * (x2 * z1 - x1 * z2)^2 
*/

void 
fmpz_ec_add(fmpz_t x, fmpz_t z, fmpz_t x1, fmpz_t z1, fmpz_t x2, fmpz_t z2,
       fmpz_t x0, fmpz_t z0, fmpz_t n)
{

    fmpz_t u, v, w;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(w);

    fmpz_sub(u, x2, z2);    /* u = (x2 - z2) */
    fmpz_add(v, x1, z1);    /* v = (x1 + z1) */

    if (fmpz_cmp(v, n) > 0)
        fmpz_sub(v, v, n);

    fmpz_mul(u, u, v);      /* u = (x2 - z2) * (x1 + z1) */
    fmpz_mod(u, u, n);

    fmpz_sub(v, x1, z1);    /* v = (x1 - z1) */
    fmpz_add(w, x2, z2);    /* w = (x2 + z2) */
    if (fmpz_cmp(w, n) > 0)
        fmpz_sub(w, w, n);

    fmpz_mul(v, v, w);      /* v = (x1 - z1) * (x2 + z2) */
    fmpz_mod(v, v, n);

    fmpz_add(w, u, v);      /* w = 2 * (x1 * x2 - z1 * z2) */
    if (fmpz_cmp(w, n) > 0)
        fmpz_sub(w, w, n);

    fmpz_sub(v, u, v);      /* v = 2 * (x2 * z1 - x1 * z2) */

    fmpz_mul(w, w, w);      /* w = 4 * (x1 * x2 - z1 * z2)^2 */
    fmpz_mod(w, w, n);

    fmpz_mul(v, v, v);      /* v = 4 * (x2 * z1 - x1 * z2)^2 */
    fmpz_mod(v, v, n);
    
    fmpz_mul(x, z0, w);     /* x = 4 * z0 * (x1 * x2 - z1 * z2)^2 */
    fmpz_mod(x, x, n);

    fmpz_mul(z, x0, v);     /* z = 4 * x0 * (x2 * z1 - x1 * z2)^2 */
    fmpz_mod(z, z, n);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(w);

}

