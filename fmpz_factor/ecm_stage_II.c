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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpn_extras.h"

/* Implementation of the stage II of ECM */

int
fmpz_factor_ecm_stage_II(mp_ptr f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                          mp_ptr n, ecm_t ecm_inf)
{
    
    mp_ptr Qx, Qz, Rx, Rz, Qdx, Qdz, a, b, g;
    mp_limb_t mmin, mmax, maxj, sz, gcdlimbs;
    int i, j, ret;
    mp_ptr arrx, arrz, Q0x2, Q0z2;

    TMP_INIT;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 

    TMP_START;
    Qx   = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Qz   = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Rx   = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Rz   = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Qdx  = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Qdz  = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Q0x2 = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    Q0z2 = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    a    = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    b    = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    g    = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    arrx = flint_malloc(((maxj >> 1) + 1) * ecm_inf->n_size * sizeof(mp_limb_t));
    arrz = flint_malloc(((maxj >> 1) + 1) * ecm_inf->n_size * sizeof(mp_limb_t));

    mpn_zero(arrx, ((maxj >> 1) + 1) * ecm_inf->n_size);
    mpn_zero(arrz, ((maxj >> 1) + 1) * ecm_inf->n_size);
    mpn_zero(Q0x2, ecm_inf->n_size);
    mpn_zero(Q0z2, ecm_inf->n_size);
    mpn_zero(g, ecm_inf->n_size);

    g[0] = ecm_inf->one[0];

    ret = 0;

    /* arr[0] = Q0 */
    flint_mpn_copyi(arrx, ecm_inf->x, ecm_inf->n_size);
    flint_mpn_copyi(arrz, ecm_inf->z, ecm_inf->n_size);

    /* Q0x2, Q0z2 = 2Q0 */
    fmpz_factor_ecm_double(Q0x2, Q0z2, arrx, arrz, n, ecm_inf);

    /* arr[1] = 3Q0 */
    fmpz_factor_ecm_add(arrx + 1 * ecm_inf->n_size, arrz + 1 * ecm_inf->n_size,
                         Q0x2, Q0z2, arrx, arrz, arrx, arrz, n, ecm_inf);

    /* For each odd j (j > 3) , compute j * Q0 [x0 :: z0] */
    /* jth stored in arr[j/2] */

    /* We are adding 2Q0 every time. Need to calculate all j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    for (j = 2; j <= (maxj >> 1); j += 1)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Difference is (j - 4)Q0 */

        fmpz_factor_ecm_add(arrx + j * ecm_inf->n_size, arrz + j * ecm_inf->n_size,
                             arrx + (j - 1) * ecm_inf->n_size, arrz + (j - 1) * ecm_inf->n_size,
                             Q0x2, Q0z2,
                             arrx + (j - 2) * ecm_inf->n_size, arrz + (j - 2) * ecm_inf->n_size,
                             n, ecm_inf);
    }

    /* Q = P * Q_0 */
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, ecm_inf->x, ecm_inf->z,
                                           P, n, ecm_inf);
    /* R = mmin * Q */
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, mmin, n, ecm_inf);

    /* Qd = (mmin - 1) * Q */
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, mmin - 1, n, ecm_inf);

    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            if (ecm_inf->prime_table[i - mmin][j] == 1)
            {
                flint_mpn_mulmod_preinvn(a, Rx, arrz + (j >> 1) * ecm_inf->n_size,
                                         ecm_inf->n_size, n, ecm_inf->ninv,
                                         ecm_inf->normbits);

                flint_mpn_mulmod_preinvn(b, Rz, arrx + (j >> 1) * ecm_inf->n_size,
                                         ecm_inf->n_size, n, ecm_inf->ninv,
                                         ecm_inf->normbits);


                fmpz_factor_ecm_submod(a, a, b, n, ecm_inf->n_size);

                flint_mpn_mulmod_preinvn(g, g, a, ecm_inf->n_size, n, ecm_inf->ninv,
                                         ecm_inf->normbits);
            }
        }

        flint_mpn_copyi(a, Rx, ecm_inf->n_size);
        flint_mpn_copyi(b, Rz, ecm_inf->n_size);

        /* R = R + Q    
           difference is stored in Qd, initially (Mmin - 1)Q */

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, ecm_inf);

        flint_mpn_copyi(Qdx, a, ecm_inf->n_size);
        flint_mpn_copyi(Qdz, b, ecm_inf->n_size);
    }

    sz = ecm_inf->n_size;
    MPN_NORM(g, sz);

    if (sz == 0)
    {
        ret = 0;
        goto cleanup;
    }

    gcdlimbs = flint_mpn_gcd_full(f, n, ecm_inf->n_size, g, sz);

    /* condition one -> gcd = n_ecm->one
       condition two -> gcd = n
       if neither is true, factor found */

    if (!(gcdlimbs == 1 && f[0] == ecm_inf->one[0]) && 
        !(gcdlimbs == ecm_inf->n_size && mpn_cmp(f, n, ecm_inf->n_size) == 0))
    {
        /* Found factor in stage II */
        ret = gcdlimbs;
        goto cleanup;
    }

    cleanup:

    TMP_END;

    flint_free(arrx);
    flint_free(arrz);

    return ret;
}
