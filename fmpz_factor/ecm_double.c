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

/* P (x : z) = 2 * P1 (x0 : z0)  */

/* 
    Coordinates of P : 

        x = (x0 + z0)^2 * (x0 - z0)^2 mod n
        z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) mod n
*/

/* a24 = (a + 2) / 4 mod n */

void
fmpz_factor_ecm_double(fmpz_t x, fmpz_t z, fmpz_t x0, fmpz_t z0, fmpz_t n, fmpz_t a24)
{
    fmpz_t u, v, w;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(w);

    fmpz_add(u, x0, z0);    /* u = x0 + z0 */
    if (fmpz_cmp(u, n) > 0)
        fmpz_sub(u, u, n);
    fmpz_mul(u, u, u);      /* u = (x0 + z0)^2 */
    fmpz_mod(u, u, n);

    fmpz_sub(v, x0, z0);    /* v = x0 - z0 */
    fmpz_mul(v, v, v);      /* v = (x0 - z0)^2 */
    fmpz_mod(v, v, n);

    fmpz_mul(x, u, v);      /* x = (x0 + z0)^2 * (x0 - z0)^2 */
    fmpz_mod(x, x, n);

    fmpz_sub(w, u, v);      /* w = 4 * x0 * z0 */
    fmpz_mul(u, w, a24);    /* u = a24 * 4 * x0 * z0 */
    fmpz_mod(u, u, n);

    fmpz_add(u, u, v);      /* u = (x0 - z0)^2 + a24 * 4 * x0 * z0 */
    if (fmpz_cmp(u, n) > 0)
        fmpz_sub(u, u, n);
    fmpz_mul(z, w, u);      /* z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) */
    fmpz_mod(z, z, n);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(w);
}


