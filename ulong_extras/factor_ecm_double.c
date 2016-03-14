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
#include "ulong_extras.h"

/* P (x : z) = 2 * P1 (x0 : z0)  */

/* 
    Coordinates of P : 

        x = (x0 + z0)^2 * (x0 - z0)^2 mod n
        z = 4 * x0 * z0 * ((x0 - z0)^2 + a24 * 4 * x0 * z0) mod n
*/

/* a24 = (a + 2) / 4 mod n */

void
n_factor_ecm_double(mp_limb_t *x, mp_limb_t *z, mp_limb_t x0, mp_limb_t z0,
                    mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w;

    if (z0 == 0)
    {
        *x = x0;
        *z = 0;
        return;
    }

    u = n_addmod(x0, z0, n);
    u = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    v = n_submod(x0, z0, n);
    v = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    *x = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    w = n_submod(u, v, n);
    u = n_mulmod_preinv(w, n_ecm_inf->a24, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    u = n_addmod(u, v, n);
    *z = n_mulmod_preinv(w, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
}
