/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Computes f such that \sigma_x(g) = f.
*/
void
unity_zp_aut(unity_zp f, const unity_zp g, ulong x)
{
    ulong i, p_pow, p_pow_preinv;
    fmpz_t coeff;
    fmpz_init(coeff);

    p_pow = n_pow(f->p, f->exp);
    p_pow_preinv = n_preinvert_limb(p_pow);

    unity_zp_set_zero(f);

    /* for i = 0, 1,..., p^k set f[i * x mod p^k] = g[i] */
    for (i = 0; i < p_pow; i++)
    {
        /* compute x * i mod p^k */
        ulong ind = n_mulmod2_preinv(x, i, p_pow, p_pow_preinv);

        /* set f[ind] = g[i] */
        fmpz_mod_poly_get_coeff_fmpz(coeff, g->poly, i, g->ctx);
        unity_zp_coeff_add_fmpz(f, ind, coeff);
    }

    _unity_zp_reduce_cyclotomic(f);
    fmpz_clear(coeff);
}

