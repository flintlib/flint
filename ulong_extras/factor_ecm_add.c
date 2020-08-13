/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"


/* P (x : z) = P1 (x1 : z1) + P2 (x2 : z2) where P0 (x0 : zo) is P - Q */

/*    Coordinates of P : 

        x = 4 * z0 * (x1 * x2 - z1 * z2)^2 mod n
        z = 4 * x0 * (x2 * z1 - x1 * z2)^2 mod n
*/

void 
n_factor_ecm_add(mp_limb_t *x, mp_limb_t *z, mp_limb_t x1, mp_limb_t z1, 
                 mp_limb_t x2, mp_limb_t z2, mp_limb_t x0, mp_limb_t z0,
                 mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w;

    if (z1 == 0)
    {
        *x = x2;
        *z = z2;
        return;
    }

    if (z2 == 0)
    {
        *x = x1;
        *z = z1;
        return;
    }

    if (z0 == 0)
    {
        n_factor_ecm_double(x, z, x1, z1, n, n_ecm_inf);
        return;
    }

    u = n_submod(x2, z2, n);        /* u = (x2 - z2) */
    v = n_addmod(x1, z1, n);        /* v = (x1 + z1) */
    
    /* u = (x2 - z2) * (x1 + z1) */
    u = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    v = n_submod(x1, z1, n);        /* v = (x1 - z1) */
    w = n_addmod(x2, z2, n);        /* w = (x2 + z2) */
    
    /* v = (x1 - z1) * (x2 + z2) */
    v = n_mulmod_preinv(v, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    w = n_addmod(u, v, n);          /* w = 2 * (x1 * x2 - z1 * z2) */
    v = n_submod(v, u, n);          /* v = 2 * (x2 * z1 - x1 * z2) */
    
    /* w = 4 * (x1 * x2 - z1 * z2)^2 */
    w = n_mulmod_preinv(w, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    
    /* v = 4 * (x2 * z1 - x1 * z2)^2 */
    v = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    
    /* x = 4 * z0 * (x1 * x2 - z1 * z2)^2 */
    *x = n_mulmod_preinv(z0, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    
    /* z = 4 * x0 * (x2 * z1 - x1 * z2)^2 */
    *z = n_mulmod_preinv(x0, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
}
