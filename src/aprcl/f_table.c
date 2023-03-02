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
    Returns a table of the function f: [1,2,...,q-2] -> [1,2,...,q-2]
    defined by $1 - g^x \equiv  gf(x) \bmod q$, 
    where $g$ is a primitive root modulo $q$.

    f_table[x - 1] = f(x).
*/
mp_ptr
aprcl_f_table(const ulong q)
{
    int i;
    ulong g, g_pow, qinv;
    mp_ptr g_table, f_table;

    g = n_primitive_root_prime(q);
    g_table = _nmod_vec_init(q);
    f_table = _nmod_vec_init(q);
    qinv = n_preinvert_limb(q);
  
    g_pow = g;
    /* g_table[g^i mod q] = i */
    for (i = 1; i < q; i++)
    {
        g_table[g_pow] = i;
        g_pow = n_mulmod2_preinv(g_pow, g, q, qinv);
    }

    g_pow = g;
    /* f_table[i] such that g^f_table[i] = 1 - g^i mod q*/
    for (i = 1; i < q; i++)
    {
        f_table[i] = g_table[n_submod(1, g_pow, q)];
        g_pow = n_mulmod2_preinv(g_pow, g, q, qinv);
    }

    _nmod_vec_clear(g_table);
    return f_table;
}
