/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Outer wrapper for ECM 
   makes calls to stage I and stage II (two) */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

/* P (x : z) = P1 (x1 : z1) + P2 (x2 : z2) where P0 (x0 : zo) is P - Q */

/*    Coordinates of P : 

        x = 4 * z0 * (x1 * x2 - z1 * z2)^2 mod n
        z = 4 * x0 * (x2 * z1 - x1 * z2)^2 mod n
*/

void 
fmpz_factor_ecm_add(mp_ptr x, mp_ptr z, mp_ptr x1, mp_ptr z1, mp_ptr x2,
                    mp_ptr z2, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf)
{

    if (flint_mpn_zero_p(z1, ecm_inf->n_size))
    {
        flint_mpn_copyi(x, x2, ecm_inf->n_size);
        flint_mpn_copyi(z, z2, ecm_inf->n_size);
        return;
    }

    if (flint_mpn_zero_p(z2, ecm_inf->n_size))
    {
        flint_mpn_copyi(x, x1, ecm_inf->n_size);
        flint_mpn_copyi(z, z1, ecm_inf->n_size);
        return;
    }

    if (flint_mpn_zero_p(z0, ecm_inf->n_size))
    {
        fmpz_factor_ecm_double(x, z, x1, z1, n, ecm_inf);
        return;
    }

    /* u = (x2 - z2) */
    fmpz_factor_ecm_submod(ecm_inf->u, x2, z2, n, ecm_inf->n_size);
    /* v = (x1 + z1) */
    fmpz_factor_ecm_addmod(ecm_inf->v, x1, z1, n, ecm_inf->n_size);

    /* u = (x2 - z2) * (x1 + z1) */
    flint_mpn_mulmod_preinvn(ecm_inf->u, ecm_inf->u, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);
    /* v = (x1 - z1) */
    fmpz_factor_ecm_submod(ecm_inf->v, x1, z1, n, ecm_inf->n_size);

    /* w = (x2 + z2) */
    fmpz_factor_ecm_addmod(ecm_inf->w, x2, z2, n, ecm_inf->n_size);
    /* v = (x1 - z1) * (x2 + z2) */
    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->v, ecm_inf->w, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);
    
    /* w = 2 * (x1 * x2 - z1 * z2) */
    fmpz_factor_ecm_addmod(ecm_inf->w, ecm_inf->u, ecm_inf->v, n, ecm_inf->n_size);
    /* v = 2 * (x2 * z1 - x1 * z2) */
    fmpz_factor_ecm_submod(ecm_inf->v, ecm_inf->v, ecm_inf->u, n, ecm_inf->n_size);

    /* w = 4 * (x1 * x2 - z1 * z2)^2 */
    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->w, ecm_inf->w, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);
    /* v = 4 * (x2 * z1 - x1 * z2)^2 */
    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->v, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    /* x = 4 * z0 * (x1 * x2 - z1 * z2)^2 */
    flint_mpn_mulmod_preinvn(x, z0, ecm_inf->w, ecm_inf->n_size, n, ecm_inf->ninv,
                             ecm_inf->normbits);
    /* z = 4 * x0 * (x2 * z1 - x1 * z2)^2 */
    flint_mpn_mulmod_preinvn(z, x0, ecm_inf->v, ecm_inf->n_size, n, ecm_inf->ninv,
                             ecm_inf->normbits);   
}
