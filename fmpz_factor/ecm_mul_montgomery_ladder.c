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

/* P (x0 : z0) <- kP using Montgomery ladder algorithm */

/* tstbit uses 0 based indexing */

void
fmpz_factor_ecm_mul_montgomery_ladder(fmpz_t x, fmpz_t z, fmpz_t x0, fmpz_t z0,
                                      fmpz_t k, fmpz_t n, ecm_t ecm_inf)
{
    mp_limb_t i;
    fmpz_t x1, z1, x2, z2;      /* Q (x1 : z1), P (x2 : z2) */

    if (fmpz_cmp_ui(k, 0) == 0)
    {
        fmpz_set_ui(x, 0);
        fmpz_set_ui(z, 0);
        return;
    }

    if (fmpz_is_one(k))
    {
        fmpz_set(x, x0);
        fmpz_set(z, z0);
        return;
    }

    fmpz_init_set(x1, x0);        /* Q <- P0 */
    fmpz_init_set(z1, z0);
    fmpz_init(x2);
    fmpz_init(z2);

    fmpz_factor_ecm_double(x2, z2, x0, z0, n, ecm_inf); /* P <- 2P0 */
    i = fmpz_sizeinbase(k, 2) - 2;

    while (1)
    {
        if (fmpz_tstbit(k, i) == 1)   /* ith bit is 1 */
        {
            fmpz_factor_ecm_add(x1, z1, x1, z1, x2, z2, x0, z0, n, ecm_inf); /* Q <- P + Q */
            fmpz_factor_ecm_double(x2, z2, x2, z2, n, ecm_inf);         /* P <- 2 * P */
        }
        else
        {   
            fmpz_factor_ecm_add(x2, z2, x1, z1, x2, z2, x0, z0, n, ecm_inf); /* P <- P + Q */
            fmpz_factor_ecm_double(x1, z1, x1, z1, n, ecm_inf);         /* Q <- 2 * Q */
        }

        if (i == 0)
            break;
        else
            i -= 1;
    }

    fmpz_set(x, x1);
    fmpz_set(z, z1);

    fmpz_clear(x1);
    fmpz_clear(z1);
    fmpz_clear(x2);
    fmpz_clear(z2);
}
