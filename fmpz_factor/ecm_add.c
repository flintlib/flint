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

/* Outer wrapper for ECM 
   makes calls to stage I and stage II (two) */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
    
/* P (x : z) = P1 (x1 : z1) + P2 (x2 : z2) where P0 (x0 : zo) is P - Q */

/*    Coordinates of P : 

        x = 4 * z0 * (x1 * x2 - z1 * z2)^2 mod n
        z = 4 * x0 * (x2 * z1 - x1 * z2)^2 mod n
*/

void 
fmpz_factor_ecm_add(fmpz_t x, fmpz_t z, fmpz_t x1, fmpz_t z1, fmpz_t x2, fmpz_t z2,
                    fmpz_t x0, fmpz_t z0, fmpz_t n, ecm_t ecm_inf)
{
    if (fmpz_cmp_ui(z1, 0) == 0)
    {
        fmpz_set(x, x2);
        fmpz_set(z, z2);
        return;
    }

    if (fmpz_cmp_ui(z2, 0) == 0)
    {

        fmpz_set(x, x1);
        fmpz_set(z, z1);
        return;
    }

    if (fmpz_cmp_ui(z0, 0) == 0)
    {
        fmpz_factor_ecm_double(x, z, x1, z1, n, ecm_inf);
        return;
    }

    fmpz_sub(ecm_inf->u, x2, z2);    /* u = (x2 - z2) */
    fmpz_add(ecm_inf->v, x1, z1);    /* v = (x1 + z1) */

    if (fmpz_cmp(ecm_inf->v, n) > 0)
        fmpz_sub(ecm_inf->v, ecm_inf->v, n);

    fmpz_mul(ecm_inf->u, ecm_inf->u, ecm_inf->v);      /* u = (x2 - z2) * (x1 + z1) */
    fmpz_mod(ecm_inf->u, ecm_inf->u, n);

    fmpz_sub(ecm_inf->v, x1, z1);    /* v = (x1 - z1) */
    fmpz_add(ecm_inf->w, x2, z2);    /* w = (x2 + z2) */
    if (fmpz_cmp(ecm_inf->w, n) > 0)
        fmpz_sub(ecm_inf->w, ecm_inf->w, n);

    fmpz_mul(ecm_inf->v, ecm_inf->v, ecm_inf->w);      /* v = (x1 - z1) * (x2 + z2) */
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);

    fmpz_add(ecm_inf->w, ecm_inf->u, ecm_inf->v);      /* w = 2 * (x1 * x2 - z1 * z2) */
    if (fmpz_cmp(ecm_inf->w, n) > 0)
        fmpz_sub(ecm_inf->w, ecm_inf->w, n);

    fmpz_sub(ecm_inf->v, ecm_inf->v, ecm_inf->u);      /* v = 2 * (x2 * z1 - x1 * z2) */

    fmpz_mul(ecm_inf->w, ecm_inf->w, ecm_inf->w);      /* w = 4 * (x1 * x2 - z1 * z2)^2 */
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);

    fmpz_mul(ecm_inf->v, ecm_inf->v, ecm_inf->v);      /* v = 4 * (x2 * z1 - x1 * z2)^2 */
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);
    
    fmpz_mul(x, z0, ecm_inf->w);     /* x = 4 * z0 * (x1 * x2 - z1 * z2)^2 */
    fmpz_mod(x, x, n);

    fmpz_mul(z, x0, ecm_inf->v);     /* z = 4 * x0 * (x2 * z1 - x1 * z2)^2 */
    fmpz_mod(z, z, n);
}