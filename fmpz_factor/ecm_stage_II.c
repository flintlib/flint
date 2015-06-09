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

/* Implementation of the stage II of ECM */

int
fmpz_factor_ecm_stage_II(fmpz_t f, fmpz_t x0, fmpz_t z0, mp_limb_t B1,
                         mp_limb_t B2, mp_limb_t D, fmpz_t a24, fmpz_t n)
{

    fmpz_t g, tim, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, mdiff, prod;
    int i, j, ret;

    mmin = (B1 + (D/2)) / D;
    mmax = ((B2 - D/2) + D - 1)/D;      /* ceil */
    maxj = (D + 1)/2; 
    mdiff = mmax - mmin + 1;

    uint8_t GCD_table[maxj + 1], prime_table[mdiff][maxj + 1];
    fmpz_t arrx[maxj + 1], arrz[maxj + 1];

    fmpz_init(Qx);
    fmpz_init(Qz);
    fmpz_init(Qdx);
    fmpz_init(Qdz);
    fmpz_init(Rx);
    fmpz_init(Rz);
    fmpz_init(a);
    fmpz_init(b);

    fmpz_init_set_ui(g, 1);

    for (j = 1; j <= maxj; j++)
    {
        fmpz_init(arrx[j]);
        fmpz_init(arrz[j]);
    }

    ret = 0;

    /* compute GCD_table */
    for (j = 1; j <= maxj; j ++)
    {
        if ((j%2) && n_gcd(j, D) == 1)
            GCD_table[j] = 1;  
        else
            GCD_table[j] = 0;
    }


    /* compute prime table */
    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            prime_table[i][j] = 0;

            /* if (i + mmin)*D + j or (i + mmin)*D - j 
               is prime, mark 1. Can be possibly prime
               only if gcd(j, D) = 1 */

            if (GCD_table[j] == 1)
            {
                prod = (i + mmin)*D + j;
                if (n_is_prime(prod))
                    prime_table[i][j] = 1;

                prod = (i + mmin)*D - j;
                if (n_is_prime(prod))
                    prime_table[i][j] = 1;
            }
        }
    }

    /* For each odd j, compute j * Q0 [x0 :: z0] */

    /* We are adding 2Q0 every time. Need to calculate all j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    /* arr[1] = Q0 */
    fmpz_set(arrx[1], x0);  
    fmpz_set(arrz[1], z0);

    /* arr[2] = 2Q0 */
    fmpz_factor_ecm_double(arrx[2], arrz[2], arrx[1], arrz[1], a24, n);

    /* arr[3] = 3Q0 */
    fmpz_factor_ecm_add(arrx[3], arrz[3], arrx[2], arrz[2], arrx[1], arrz[1], 
                        arrx[1], arrz[1], a24, n);


    for (j = 5; j <= maxj; j += 2)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Differnce is (j - 4)Q0 */

        fmpz_factor_ecm_add(arrx[j], arrz[j], arrx[j - 2], arrz[j - 2], 
                            arrx[2], arrz[2], arrx[j - 4], arrz[j - 4], 
                            a24, n);
    }


    /* Q = D * Q_0 */
    fmpz_set_ui(tim, D);
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, x0, z0, tim, a24, n);
    
    /* R = mmin * Q */
    fmpz_set_ui(tim, mmin);
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, tim, a24, n);

    /* Qd = (mmin - 1) * Q */
    fmpz_set_ui(tim, mmin - 1);
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, tim, a24, n);
                
    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
        for (j = 1; j <= maxj; j+=2)
        {
            if ((GCD_table[j] == 1) && (prime_table[i - mmin][j] == 1))
            {
                fmpz_mul(a, Rx, arrz[j]);
                fmpz_mod(a, a, n);

                fmpz_mul(b, Rz, arrx[j]);
                fmpz_mod(b, b, n);

                fmpz_sub(a, a, b);
                fmpz_mul(g, g, a);
                fmpz_mod(g, g, n);
            }
        }

        fmpz_set(a, Rx);
        fmpz_set(b, Rz);

        /* R = R + Q    
           difference is stored in Qd, initially (Mmin - 1)Q */

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, a24, n);

        fmpz_set(Qdx, a);
        fmpz_set(Qdz, b);
    }

    fmpz_gcd(f, g, n);

    if (!fmpz_is_one(f) && fmpz_cmp(f, n))
        ret = 1;

    fmpz_clear(Qx);
    fmpz_clear(Qz);
    fmpz_clear(Qdx);
    fmpz_clear(Qdz);
    fmpz_clear(Rx);
    fmpz_clear(Rz);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);

    for (j = 1; j <= maxj; j++)
    {
        fmpz_clear(arrx[j]);
        fmpz_clear(arrz[j]);
    }

    return ret;
}
