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
    mp_size_t invlimbs;
    mp_ptr temp;
    int ret;
    fmpz_t a, b, ff, inv;
    __mpz_struct *aa, *bb;

    ret = 0;
    temp = flint_malloc(ecm_inf->n_size * sizeof(mp_limb_t));

    mpn_zero(temp, ecm_inf->n_size);
    mpn_copyi(ecm_inf->u, sig, ecm_inf->n_size);

    temp[0] = UWORD(4);
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

    fmpz_init(ff);
    fmpz_init(inv);
    fmpz_init2(a, ecm_inf->n_size);
    fmpz_init2(b, ecm_inf->n_size);

    aa = _fmpz_promote(a);
    bb = _fmpz_promote(b);

    mpn_copyi(aa->_mp_d, ecm_inf->v, sz);
    mpn_copyi(bb->_mp_d, n, ecm_inf->n_size);

    aa->_mp_size = sz;
    bb->_mp_size = ecm_inf->n_size;

    mpn_rshift(aa->_mp_d, aa->_mp_d, aa->_mp_size, ecm_inf->normbits);
    mpn_rshift(bb->_mp_d, bb->_mp_d, bb->_mp_size, ecm_inf->normbits);

    fmpz_gcdinv(ff, inv, a, b);

    if (fmpz_is_one(ff) == 0)
    {
        aa = _fmpz_promote(ff);
        MPN_NORM(aa->_mp_d, aa->_mp_size);

        mpn_zero(f, ecm_inf->n_size);
        mpn_copyi(f, aa->_mp_d, aa->_mp_size);

        ret = aa->_mp_size;
        goto cleanup;
    }

    fmpz_mod(inv, inv, b);

    bb = _fmpz_promote(inv);
    mpn_zero(ecm_inf->u, ecm_inf->n_size);

    cy = mpn_lshift(ecm_inf->u, bb->_mp_d, bb->_mp_size, ecm_inf->normbits);

    if (cy)
        ecm_inf->u[bb->_mp_size] = cy;

    invlimbs = bb->_mp_size;

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

    flint_free(temp);

    fmpz_clear(ff);
    fmpz_clear(inv);
    fmpz_clear(a);
    fmpz_clear(b);
    return ret;
}
