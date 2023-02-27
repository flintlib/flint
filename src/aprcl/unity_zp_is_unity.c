/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

slong
unity_zp_is_unity(unity_zp f)
{
    ulong result;
    ulong i, p_pow;
    unity_zp unity;

    p_pow = n_pow(f->p, f->exp);
    unity_zp_init(unity, f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));

    /* if the power was not found returns -1 */
    result = -1;
    for (i = 0; i < p_pow; i++)
    {
        /* set unity = \zeta_{p^k}^i */
        unity_zp_set_zero(unity);
        unity_zp_coeff_set_ui(unity, i, 1);

        /* check if f = zeta_{p^k}^i */
        if (unity_zp_equal(unity, f) == 1)
        {
            /* if so, returns \zeta_{p^k} power */
            result = i;
            break;
        }
    }

    unity_zp_clear(unity);
    return result;
}

