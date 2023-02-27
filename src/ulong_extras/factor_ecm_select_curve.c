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
#include "ulong_extras.h"

int
n_factor_ecm_select_curve(mp_limb_t *f, mp_limb_t sig, mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t u, v, w, t, hi, lo;
    mp_ptr a;
    TMP_INIT;

    TMP_START;
    a = TMP_ALLOC(2 * sizeof(mp_limb_t));

    u = sig;

    /* v = sig * 4 */
    v = n_mulmod_preinv(u, UWORD(4) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);

    /* w = sig ^ 2 */
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    u = w - (UWORD(5) << n_ecm_inf->normbits);  /* u = sig^2 - 5 */

    /* w = u * u */
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* x = u * u * u */
    n_ecm_inf->x = n_mulmod_preinv(w, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    /* w = v * v */
    w = n_mulmod_preinv(v, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits); 

    /* z = v * v * v */
    n_ecm_inf->z = n_mulmod_preinv(w, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    w = n_mulmod_preinv(n_ecm_inf->x, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    t = n_mulmod_preinv(w, UWORD(4) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);
    w = n_mulmod_preinv(u, UWORD(3) << n_ecm_inf->normbits, n, n_ecm_inf->ninv,
                        n_ecm_inf->normbits);
    
    u = n_submod(v, u, n);
    v = n_addmod(v, w, n);
    w = n_mulmod_preinv(u, u, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    u = n_mulmod_preinv(u, w, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    
    n_ecm_inf->a24 = n_mulmod_preinv(u, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    v = n_mulmod_preinv(t, n_ecm_inf->z, n, n_ecm_inf->ninv, n_ecm_inf->normbits);

    *f = n_gcdinv(&u, v, n);

    if (*f == n)
        return 0;
    else if (*f > n_ecm_inf->one)
        return 1;

    a[1] = UWORD(0);
    a[0] = u;
    mpn_lshift(a, a, 2, n_ecm_inf->normbits);       /* shifting back right */
    hi = a[1];
    lo = a[0];
    u = n_ll_mod_preinv(hi, lo, n, n_ecm_inf->ninv);

    v = n_mulmod_preinv(u, t, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    n_ecm_inf->x = n_mulmod_preinv(n_ecm_inf->x, v, n, n_ecm_inf->ninv, 
                                   n_ecm_inf->normbits);

    v = n_mulmod_preinv(u, n_ecm_inf->z, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    w = n_mulmod_preinv(n_ecm_inf->a24, v, n, n_ecm_inf->ninv, n_ecm_inf->normbits);
    n_ecm_inf->a24 = n_addmod(w, UWORD(2) << n_ecm_inf->normbits, n);
    n_ecm_inf->a24 >>= 2;
    n_ecm_inf->a24 >>= n_ecm_inf->normbits;
    n_ecm_inf->a24 <<= n_ecm_inf->normbits;
    n_ecm_inf->z = n_ecm_inf->one;

    TMP_END;

    return 0;
}
