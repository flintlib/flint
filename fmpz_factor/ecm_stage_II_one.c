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
#include "fmpz_vec.h"

/* Implementation of the stage I of ECM */

int
fmpz_factor_ecm_stage_II_one(fmpz_t f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                             fmpz_t n, ecm_t ecm_inf)
{

    fmpz_t g, tim, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj;
    int i, j, ret;
    fmpz * arrx, * arrz;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 

    fmpz_init(Qx);
    fmpz_init(Qz);
    fmpz_init(Qdx);
    fmpz_init(Qdz);
    fmpz_init(Rx);
    fmpz_init(Rz);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(tim);

    fmpz_init_set_ui(g, 1);

    arrx = _fmpz_vec_init(maxj + 1);
    arrz = _fmpz_vec_init(maxj + 1);

    ret = 0;

    /* arr[1] = Q0 */
    fmpz_set(arrx + 1, ecm_inf->x);  
    fmpz_set(arrz + 1, ecm_inf->z);

    /* arr[2] = 2Q0 */
    fmpz_factor_ecm_double(arrx + 2, arrz + 2, arrx + 1, arrz + 1, n, ecm_inf);

    /* arr[3] = 3Q0 */
    fmpz_factor_ecm_add(arrx + 3, arrz + 3, arrx + 2, arrz + 2, arrx + 1, arrz + 1, 
                        arrx + 1, arrz + 1, n, ecm_inf);

    /* For each odd j (j > 3) , compute j * Q0 [x0 :: z0] */

    /* We are adding 2Q0 every time. Need to calculate all j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    for (j = 5; j <= maxj; j += 2)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Differnce is (j - 4)Q0 */

        fmpz_factor_ecm_add(arrx + j, arrz + j, arrx + j - 2, arrz + j - 2, 
                            arrx + 2, arrz + 2, arrx + j - 4, arrz + j - 4, 
                            n, ecm_inf);
    }

    /* Q = D * Q_0 */
    fmpz_set_ui(tim, P);
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, ecm_inf->x, ecm_inf->z, tim, n, ecm_inf);
    
    /* R = mmin * Q */
    fmpz_set_ui(tim, mmin);
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, tim, n, ecm_inf);

    /* Qd = (mmin - 1) * Q */
    fmpz_set_ui(tim, mmin - 1);
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, tim, n, ecm_inf);
                
    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
        for (j = 1; j <= maxj; j+=2)
        {
            if (ecm_inf->prime_table[i - mmin][j] == 1)
            {
                fmpz_mul(a, Rx, arrz + j);
                fmpz_mod(a, a, n);

                fmpz_mul(b, Rz, arrx + j);
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

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, ecm_inf);

        fmpz_set(Qdx, a);
        fmpz_set(Qdz, b);
    }


    fmpz_gcd(f, g, n);

    if (!fmpz_is_one(f) && fmpz_cmp(f, n))
        ret = 1;

    fmpz_clear(tim);
    fmpz_clear(Qx);
    fmpz_clear(Qz);
    fmpz_clear(Qdx);
    fmpz_clear(Qdz);
    fmpz_clear(Rx);
    fmpz_clear(Rz);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);

    _fmpz_vec_clear(arrx, maxj + 1);
    _fmpz_vec_clear(arrz, maxj + 1);

    return ret;
}
