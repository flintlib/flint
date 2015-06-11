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
fmpz_factor_ecm_double(fmpz_t x, fmpz_t z, fmpz_t x0, fmpz_t z0, fmpz_t n, ecm_t ecm_inf)
{
    if (fmpz_cmp_ui(z0, 0) == 0)
    {
        fmpz_set(x, x0);
        fmpz_set(z, z0);
        return;
    }

    fmpz_add(ecm_inf->u, x0, z0);    /* u = x0 + z0 */
    if (fmpz_cmp(ecm_inf->u, n) > 0)
        fmpz_sub(ecm_inf->u, ecm_inf->u, n);
    fmpz_mul(ecm_inf->u, ecm_inf->u, ecm_inf->u);      /* u = (x0 + z0)^2 */
    fmpz_mod(ecm_inf->u, ecm_inf->u, n);

    fmpz_sub(ecm_inf->v, x0, z0);    /* v = x0 - z0 */
    fmpz_mul(ecm_inf->v, ecm_inf->v, ecm_inf->v);      /* v = (x0 - z0)^2 */
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);

    fmpz_mul(x, ecm_inf->u, ecm_inf->v);      /* x = (x0 + z0)^2 * (x0 - z0)^2 */
    fmpz_mod(x, x, n);

    fmpz_sub(ecm_inf->w, ecm_inf->u, ecm_inf->v);      /* w = 4 * x0 * z0 */
    fmpz_mul(ecm_inf->u, ecm_inf->w, ecm_inf->a24);    /* u = a24 * 4 * x0 * z0 */
    fmpz_mod(ecm_inf->u, ecm_inf->u, n);

    fmpz_add(ecm_inf->u, ecm_inf->u, ecm_inf->v);      /* u = (x0 - z0)^2 + a24 * 4 * x0 * z0 */
    if (fmpz_cmp(ecm_inf->u, n) > 0)
        fmpz_sub(ecm_inf->u, ecm_inf->u, n);
    fmpz_mul(z, ecm_inf->w, ecm_inf->u);      /* z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) */
    fmpz_mod(z, z, n);
}
