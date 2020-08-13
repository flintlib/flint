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
#include "nmod_vec.h"
#include "ulong_extras.h"

int
n_factor_ecm_stage_II(mp_limb_t *f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                      mp_limb_t n, n_ecm_t n_ecm_inf)
{
    
    mp_limb_t g, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, Q0x2, Q0z2;
    int i, j, ret;
    mp_ptr arrx, arrz;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 

    g = n_ecm_inf->one;

    arrx = _nmod_vec_init((maxj >> 1) + 1);
    arrz = _nmod_vec_init((maxj >> 1) + 1);

    ret = 0;

    /* arr[0] = Q0 */
    arrx[0] = n_ecm_inf->x;
    arrz[0] = n_ecm_inf->z;

    /* Q0x2, Q0z2 = 2Q0 */
    n_factor_ecm_double(&Q0x2, &Q0z2, arrx[0], arrz[0], n, n_ecm_inf);

    /* arr[1] = 3Q0 */
    n_factor_ecm_add(arrx + 1, arrz + 1, Q0x2, Q0z2, arrx[0], arrz[0],
                     arrx[0], arrz[0], n, n_ecm_inf);

    /* For each odd j (j > 3) , compute j * Q0 [x0 :: z0] */
    /* jth stored in arr[j/2] */

    /* We are adding 2Q0 every time. Need to calculate all j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    for (j = 2; j <= (maxj >> 1); j += 1)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Differnce is (j - 4)Q0 */

        n_factor_ecm_add(arrx + j, arrz + j, arrx[j - 1], arrz[j - 1], Q0x2,
                         Q0z2, arrx[j - 2], arrz[j - 2], n, n_ecm_inf);
    }

    /* Q = D * Q_0 */
    n_factor_ecm_mul_montgomery_ladder(&Qx, &Qz, n_ecm_inf->x, n_ecm_inf->z, P, n, n_ecm_inf);

    /* R = mmin * Q */
    n_factor_ecm_mul_montgomery_ladder(&Rx, &Rz, Qx, Qz, mmin, n, n_ecm_inf);

    /* Qd = (mmin - 1) * Q */
    n_factor_ecm_mul_montgomery_ladder(&Qdx, &Qdz, Qx, Qz, mmin - 1, n, n_ecm_inf);

    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            if (n_ecm_inf->prime_table[i - mmin][j] == 1)
            {
                a = n_mulmod_preinv(Rx, arrz[j >> 1], n, n_ecm_inf->ninv, n_ecm_inf->normbits);
                b = n_mulmod_preinv(Rz, arrx[j >> 1], n, n_ecm_inf->ninv, n_ecm_inf->normbits);
                a = n_submod(a, b, n);
                g = n_mulmod_preinv(g, a, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
            }
        }

        a = Rx;
        b = Rz;

        /* R = R + Q    
           difference is stored in Qd, initially (Mmin - 1)Q */

        n_factor_ecm_add(&Rx, &Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, n_ecm_inf);

        Qdx = a;
        Qdz = b;
    }

    *f = n_gcd(g, n);

    if ((*f > n_ecm_inf->one) && (*f < n))
    {
        /* Found factor in stage I */
        ret = 1;
    }

    _nmod_vec_clear(arrx);
    _nmod_vec_clear(arrz);

    return ret;
}
