/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

/* P (x0 : z0) <- kP using Montgomery ladder algorithm */

void
n_factor_ecm_mul_montgomery_ladder(ulong_ptr x, ulong_ptr z, ulong x0, ulong z0,
                                   ulong k, ulong n, n_ecm_t n_ecm_inf)
{
    ulong x1, z1, x2, z2, len;      /* Q (x1 : z1), P (x2 : z2) */

    if (k == 0)
    {
        *x = 0;
        *z = 0;
        return;
    }

    if (k == 1)
    {
        *x = x0;
        *z = z0;
        return;
    }

    x1 = x0;    /* Q <- P0 */
    z1 = z0;
    x2 = 0;
    z2 = 0;

    n_factor_ecm_double(&x2, &z2, x0, z0, n, n_ecm_inf);    /* P <- 2P0 */

    len = n_sizeinbase(k, 2) - 2;

    while (1)
    {
        if (((UWORD(1) << len) & k) != 0)       /* ith bit is 1 */
        {
            /* Q <- P + Q */
            n_factor_ecm_add(&x1, &z1, x1, z1, x2, z2, x0, z0, n, n_ecm_inf);   
            /* P <- 2 * P */
            n_factor_ecm_double(&x2, &z2, x2, z2, n, n_ecm_inf);                
        }
        else
        {   
            /* P <- P + Q */
            n_factor_ecm_add(&x2, &z2, x1, z1, x2, z2, x0, z0, n, n_ecm_inf);   
            /* Q <- 2 * Q */        
            n_factor_ecm_double(&x1, &z1, x1, z1, n, n_ecm_inf);                
        }


        if (len == 0)
            break;
        else
            len -= 1;
    }

    *x = x1;
    *z = z1;
}
