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
    Computes gauss sum for character \chi corresponding (q, p).
*/
void
unity_zpq_gauss_sum(unity_zpq f, ulong q, ulong p)
{
    slong i, qinv, qpow, ppow, g;

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(q);
    qpow = 1;
    ppow = 0;

    for (i = 1; i < q; i++)
    {
        qpow = n_mulmod2_preinv(qpow, g, q, qinv);
        ppow = n_addmod(ppow, 1, p);
        unity_zpq_coeff_add_ui(f, qpow, ppow, 1);
    }
}

