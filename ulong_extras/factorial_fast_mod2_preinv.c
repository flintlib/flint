/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_limb_t
n_factorial_fast_mod2_preinv(ulong n, mp_limb_t p, mp_limb_t pinv)
{
    slong i, m;
    nmod_t mod;
    mp_ptr t, u, v;
    mp_limb_t r, s;

    if (p == UWORD(1) || n >= p)
        return UWORD(0);

    if (n <= 1)
        return UWORD(1);

    nmod_init(&mod, p);

    m = n_sqrt(n);

    t = _nmod_vec_init(m + 1);
    u = _nmod_vec_init(m + 1);
    v = _nmod_vec_init(m + 1);

    t[0] = UWORD(0);
    for (i = 1; i < m; i++)
        t[i] = n_submod(t[i-1], UWORD(1), p);

    _nmod_poly_product_roots_nmod_vec(u, t, m, mod);

    for (i = 0; i < m; i++)
        t[i] = n_mod2_preinv(i * m + 1, p, pinv);

    _nmod_poly_evaluate_nmod_vec_fast(v, u, m + 1, t, m, mod);

    r = 1;
    for (i = 0; i < m; i++)
        r = n_mulmod2_preinv(r, v[i], mod.n, mod.ninv);

    for (s = m * m + 1; s <= n; s++)
        r = n_mulmod2_preinv(r, s, mod.n, mod.ninv);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
    _nmod_vec_clear(v);

    return r;
}

