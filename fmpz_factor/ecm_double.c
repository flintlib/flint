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
#include "mpn_extras.h"

/* P (x : z) = 2 * P1 (x0 : z0)  */

/* 
    Coordinates of P : 

        x = (x0 + z0)^2 * (x0 - z0)^2 mod n
        z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) mod n
*/

void
fmpz_factor_ecm_double(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0,
                       mp_ptr n, ecm_t ecm_inf)
{
    if (flint_mpn_zero_p(z0, ecm_inf->n_size))
    {
        mpn_copyi(x, x0, ecm_inf->n_size);
        mpn_zero(z, ecm_inf->n_size);
        return;
    }

    /* u = x0 + z0 */
    fmpz_factor_ecm_addmod(ecm_inf->u, x0, z0, n, ecm_inf->n_size);
    
    /* u = (x0 + z0)^2 */
    flint_mpn_mulmod_preinvn(ecm_inf->u, ecm_inf->u, ecm_inf->u, ecm_inf->n_size, n,
                             ecm_inf->ninv, ecm_inf->normbits);

    /* v = x0 - z0 */
    fmpz_factor_ecm_submod(ecm_inf->v, x0, z0, n, ecm_inf->n_size);
        
    /* v = (x0 - z0)^2 */
    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->v, ecm_inf->v, ecm_inf->n_size, n,
                             ecm_inf->ninv, ecm_inf->normbits);

    /* x = (x0 + z0)^2 * (x0 - z0)^2 */
    flint_mpn_mulmod_preinvn(x, ecm_inf->u, ecm_inf->v, ecm_inf->n_size, n,
                             ecm_inf->ninv, ecm_inf->normbits);
    /* w = 4 * x0 * z0 */
    fmpz_factor_ecm_submod(ecm_inf->w, ecm_inf->u, ecm_inf->v, n, ecm_inf->n_size);

    /* u = a24 * 4 * x0 * z0 */
    flint_mpn_mulmod_preinvn(ecm_inf->u, ecm_inf->w, ecm_inf->a24, ecm_inf->n_size, n,
                             ecm_inf->ninv, ecm_inf->normbits);
    
    /* u = (x0 - z0)^2 + a24 * 4 * x0 * z0 */
    fmpz_factor_ecm_addmod(ecm_inf->u, ecm_inf->u, ecm_inf->v, n, ecm_inf->n_size);

    /* z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) */
    flint_mpn_mulmod_preinvn(z, ecm_inf->w, ecm_inf->u, ecm_inf->n_size, n,
                             ecm_inf->ninv, ecm_inf->normbits);
}
