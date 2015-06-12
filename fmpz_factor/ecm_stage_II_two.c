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

/* Implementation of the stage I of ECM */

int
fmpz_factor_ecm_stage_II_two(fmpz_t f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                             fmpz_t n, ecm_t ecm_inf)
{
    fmpz_t g, tim, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, mdiff, prod, Psq;
    int i, j, ret;
    fmpz *arrx, *arrz;

    fmpz_init(Qx);
    fmpz_init(Qz);
    fmpz_init(Qdx);
    fmpz_init(Qdz);
    fmpz_init(Rx);
    fmpz_init(Rz);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init_set_ui(g, 1);

    Psq = P*P;
    mmin = (B1) / (Psq);
    mmax = (B2) / (Psq);
    mdiff = mmax - mmin + 1;

    ret = 0;
    j = 0;

    fmpz_t jx[Psq], jz[Psq];

    fmpz_init(jx[2]);
    fmpz_init(jz[2]);

    for (i = 1; i < P; i += 2)
    {
        fmpz_init(jx[i]);
        fmpz_init(jz[i]);
    }

    /* For each odd j (j < P) , compute j * Q0 [x0 :: z0]  */

    /* We are adding 2Q0 every time. Need to calculate all odd j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    /* arr[1] = Q0 */
    fmpz_set(jx[1], ecm_inf->x);  
    fmpz_set(jz[1], ecm_inf->z);

    /* arr[2] = 2Q0 */
    fmpz_factor_ecm_double(jx[2], jz[2], jx[1], jz[1], n, ecm_inf);

    /* arr[3] = 3Q0 */
    fmpz_factor_ecm_add(jx[3], jz[3], jx[2], jz[2], jx[1], jz[1], 
                        jx[1], jz[1], n, ecm_inf);

    for (j = 5; j < P; j += 2)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Differnce is (j - 4)Q0 */

        fmpz_factor_ecm_add(jx[j], jz[j], jx[j - 2], jz[j - 2], 
                            jx[2], jz[2], jx[j - 4], jz[j - 4], 
                            n, ecm_inf);     
    }

    /* arr[P] = double(arr[P/2]) */
    /* P is of form 2 * odd, hence P/2 will always exist in array (is odd)*/
    
    fmpz_init(jx[P]);
    fmpz_init(jz[P]);
    fmpz_factor_ecm_double(jx[P], jz[P], jx[P/2], jz[P/2], n, ecm_inf);

    /* computing P + i (i < P, gcd(i, P) = 1) from P, i, (P - i) */

    for (i = P + 1; i < 2*P; i += 2)
    {
        int val = i - P;
        if (ecm_inf->GCD_table[val] == 1)
        {
            fmpz_init(jx[i]);
            fmpz_init(jz[i]);
            fmpz_factor_ecm_add(jx[i], jz[i], jx[P], jz[P], 
                                jx[val], jz[val], jx[P - val], jz[P - val],
                                n, ecm_inf);
        }
    }

    /* computing all P + i, 2P + i, ... , (P - 1)P + i 
       s.t. (i < P) & gcd(i, P) = 1 */

    /* computing nP + i (n > 1) from (n - 1)P, P and (n - 2)P */

    for (i = 1; i < P; i += 2)
    {
        if (ecm_inf->GCD_table[i] == 1)
        {
            for (j = 2; j < P; j += 1)
            {
                fmpz_init(jx[j*P + i]);
                fmpz_init(jz[j*P + i]);            
                fmpz_factor_ecm_add(jx[j*P + i], jz[j*P + i], jx[(j-1)*P + i],
                                    jz[(j-1)*P + i], jx[P], jz[P], 
                                    jx[(j-2)*P + i], jz[(j-2)*P + i], n, 
                                    ecm_inf);
            }   
        }
    }

    /* Q = P^2 * Q0 , point to be added (giant step difference) */

    fmpz_set_ui(tim, Psq);
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, ecm_inf->x, ecm_inf->z, tim,
                                          n, ecm_inf);
    
    /* R = mmin * Q , initial point (giant step) */

    fmpz_set_ui(tim, mmin);
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, tim, n, ecm_inf);

    /* Qd = (mmin - 1) * Q , difference between current gaint step 
       and giant step difference. 
       Difference between point and difference, *inception moment* */

    fmpz_set_ui(tim, mmin - 1);
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, tim, n, ecm_inf);
           
    /* checking whether prime formed by (i*Q + j) is identity */
    for (i = mmin; i <= mmax; i += 1)
    {
        for (j = 1; j <= Psq; j+=2)
        {
            if (ecm_inf->prime_table[i - mmin][j] == 1)
            {
                fmpz_mul(a, Rx, jz[j]);
                fmpz_mod(a, a, n);

                fmpz_mul(b, Rz, jx[j]);
                fmpz_mod(b, b, n);

                fmpz_sub(a, a, b);
                fmpz_mul(g, g, a);
                fmpz_mod(g, g, n);
            }
        }

        fmpz_set(a, Rx);
        fmpz_set(b, Rz);

        /* R = R + Q    
           difference is stored in Qd, initially (mmin - 1)Q */

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, ecm_inf);

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

    return ret;
}
