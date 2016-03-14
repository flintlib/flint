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
#include "mpn_extras.h"

/* Select Montgomery Elliptic Curve given a sigma
   (Suyama's parameterization) 
   Returns 1 in case factor is found while selecting
   the curev. */

/* Also selects initial point Q0 [x0 :: z0]  (z0 = 1) */

int
fmpz_factor_ecm_select_curve(mp_ptr f, mp_ptr sig, mp_ptr n, ecm_t ecm_inf)
{
    mp_size_t sz, cy;
    mp_size_t invlimbs, gcdlimbs;
    mp_ptr temp, tempv, tempn, tempi, tempf;
    int ret;

    TMP_INIT;

    TMP_START;
    temp = TMP_ALLOC(ecm_inf->n_size * sizeof(mp_limb_t));
    tempv = TMP_ALLOC((ecm_inf->n_size) * sizeof(mp_limb_t));
    tempn = TMP_ALLOC((ecm_inf->n_size) * sizeof(mp_limb_t));
    tempi = TMP_ALLOC((ecm_inf->n_size + 1) * sizeof(mp_limb_t));
    tempf = TMP_ALLOC((ecm_inf->n_size + 1) * sizeof(mp_limb_t));

    mpn_zero(tempn, ecm_inf->n_size);
    mpn_zero(tempv, ecm_inf->n_size);
    mpn_zero(temp, ecm_inf->n_size);
    mpn_copyi(ecm_inf->u, sig, ecm_inf->n_size);

    temp[0] = UWORD(4);
    ret = 0;

    mpn_lshift(temp, temp, ecm_inf->n_size, ecm_inf->normbits);   /* temp = (4 << norm) */

    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->u, temp, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->u, ecm_inf->u, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    mpn_add_n(temp, temp, ecm_inf->one, ecm_inf->n_size); /* temp = (5 << norm) */

    fmpz_factor_ecm_submod(ecm_inf->u, ecm_inf->w, temp, n, ecm_inf->n_size);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->u, ecm_inf->u, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->x, ecm_inf->w, ecm_inf->u, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->v, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->z, ecm_inf->w, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->x, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    mpn_sub_n(temp, temp, ecm_inf->one, ecm_inf->n_size); /* temp = (4 << norm) */
    
    flint_mpn_mulmod_preinvn(ecm_inf->t, ecm_inf->w, temp, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    mpn_sub_n(temp, temp, ecm_inf->one, ecm_inf->n_size); /* temp = (3 << norm) */

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->u, temp, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    fmpz_factor_ecm_submod(ecm_inf->u, ecm_inf->v, ecm_inf->u, n, ecm_inf->n_size);

    fmpz_factor_ecm_addmod(ecm_inf->v, ecm_inf->v, ecm_inf->w, n, ecm_inf->n_size);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->u, ecm_inf->u, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->u, ecm_inf->u, ecm_inf->w, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->a24, ecm_inf->u, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->t, ecm_inf->z, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    sz = ecm_inf->n_size;
    MPN_NORM(ecm_inf->v, sz);  /* sz has number of limbs of v */

    if (sz == 0)        /* v = 0, gcd(0, n) = n. Hence inverse will not exist. */
    {                   /* No point in further computation with curve */
        ret = (-1);
        goto cleanup;
    }

    mpn_copyi(tempv, ecm_inf->v, sz);
    mpn_copyi(tempn, n, ecm_inf->n_size);

    gcdlimbs = mpn_gcdext(tempf, tempi, &invlimbs, tempv, sz, tempn, ecm_inf->n_size);

    if ((((gcdlimbs == 1) && tempf[0] == ecm_inf->one[0]) || 
        ((gcdlimbs == ecm_inf->n_size) && mpn_cmp(tempf, n, ecm_inf->n_size) == 0)) == 0)
    {
        /* Found factor */
        ret = gcdlimbs;
        goto cleanup;
    }

    if (invlimbs < 0) /* negative inverse, add n */
    {
        invlimbs *= -1;

        cy = mpn_lshift(tempi, tempi, invlimbs, ecm_inf->normbits);
        if (cy)
            tempi[invlimbs] = cy;

        mpn_sub_n(tempi, n, tempi, ecm_inf->n_size);\
    }
    else
    {
        cy = mpn_lshift(tempi, tempi, invlimbs, ecm_inf->normbits);
        if (cy)
            tempi[invlimbs] = cy;
    }

    MPN_NORM(tempi, invlimbs);
    mpn_zero(ecm_inf->u, ecm_inf->n_size);
    mpn_copyi(ecm_inf->u, tempi, invlimbs);

    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->u, ecm_inf->t, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->x, ecm_inf->x, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->v, ecm_inf->u, ecm_inf->z, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    flint_mpn_mulmod_preinvn(ecm_inf->w, ecm_inf->a24, ecm_inf->v, ecm_inf->n_size,
                             n, ecm_inf->ninv, ecm_inf->normbits);

    mpn_zero(temp, sz);
    temp[0] = UWORD(2);
    mpn_lshift(temp, temp, ecm_inf->n_size, ecm_inf->normbits);

    fmpz_factor_ecm_addmod(ecm_inf->a24, ecm_inf->w, temp, n, ecm_inf->n_size);
    mpn_rshift(ecm_inf->a24, ecm_inf->a24, ecm_inf->n_size, 2 + ecm_inf->normbits);
    mpn_lshift(ecm_inf->a24, ecm_inf->a24, ecm_inf->n_size, ecm_inf->normbits);

    mpn_copyi(ecm_inf->z, ecm_inf->one, ecm_inf->n_size);

    cleanup:

    TMP_END;
    
    return ret;
}
