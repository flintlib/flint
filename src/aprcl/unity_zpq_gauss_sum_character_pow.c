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
    Computes gauss sum for character \chi^n corresponding (q, p).
*/
void
unity_zpq_gauss_sum_character_pow(unity_zpq f, ulong q, ulong p, ulong pow)
{
    slong i, qinv, pinv, qpow, ppow, g;

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(q);
    pinv = n_preinvert_limb(p);
    qpow = 1;

    for (i = 1; i < q; i++)
    {
        qpow = n_mulmod2_preinv(qpow, g, q, qinv);
        ppow = n_mulmod2_preinv(i, pow, p, pinv);
        unity_zpq_coeff_add_ui(f, qpow, ppow, 1);
    }
}

/*
    Computes gauss sum for character \chi^n corresponding (q, p).
*/
void
unity_zpq_gauss_sum_sigma_pow(unity_zpq f, ulong q, ulong p)
{
    ulong n;

    n = fmpz_fdiv_ui(fmpz_mod_ctx_modulus(f->ctx), p);
    unity_zpq_gauss_sum_character_pow(f, q, p, n);
}

