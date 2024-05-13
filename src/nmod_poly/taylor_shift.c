/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_taylor_shift(nn_ptr poly, ulong c, slong len, nmod_t mod)
{
    if (len < 100 || (ulong) len > mod.n)
        _nmod_poly_taylor_shift_horner(poly, c, len, mod);
    else if ((c == 1 || c == mod.n - 1) && len < 1000)
        _nmod_poly_taylor_shift_horner(poly, c, len, mod);
    else
        _nmod_poly_taylor_shift_convolution(poly, c, len, mod);
}

void
nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, ulong c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift(g->coeffs, c, g->length, g->mod);
}

void
_nmod_poly_taylor_shift_convolution(nn_ptr p, ulong c, slong len, nmod_t mod)
{
    slong i, n = len - 1;
    ulong f, d;
    nn_ptr t, u;

    if (c == 0 || len <= 1)
        return;

    t = _nmod_vec_init(len);
    u = _nmod_vec_init(len);

    f = 1;
    for (i = 2; i <= n; i++)
    {
        f = n_mulmod2_preinv(f, i, mod.n, mod.ninv);
        p[i] = n_mulmod2_preinv(p[i], f, mod.n, mod.ninv);
    }

    _nmod_poly_reverse(p, p, len, len);

    t[n] = 1;
    for (i = n; i > 0; i--)
        t[i - 1] = n_mulmod2_preinv(t[i], i, mod.n, mod.ninv);

    if (c == mod.n - 1)
    {
        for (i = 1; i <= n; i += 2)
            t[i] = nmod_neg(t[i], mod);
    }
    else if (c != 1)
    {
        d = c;

        for (i = 1; i <= n; i++)
        {
            t[i] = n_mulmod2_preinv(t[i], d, mod.n, mod.ninv);
            d = n_mulmod2_preinv(d, c, mod.n, mod.ninv);
        }
    }

    _nmod_poly_mullow(u, p, len, t, len, len, mod);

    f = n_mulmod2_preinv(f, f, mod.n, mod.ninv);
    f = n_invmod(f, mod.n);

    for (i = n; i >= 0; i--)
    {
        p[i] = n_mulmod2_preinv(u[n - i], f, mod.n, mod.ninv);
        f = n_mulmod2_preinv(f, (i == 0) ? 1 : i, mod.n, mod.ninv);
    }

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_taylor_shift_convolution(nmod_poly_t g, const nmod_poly_t f,
    ulong c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift_convolution(g->coeffs, c, g->length, g->mod);
}

void
_nmod_poly_taylor_shift_horner(nn_ptr poly, ulong c, slong n, nmod_t mod)
{
    slong i, j;

    if (c == 1)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                poly[j] = nmod_add(poly[j], poly[j + 1], mod);
    }
    else if (c == mod.n - 1)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                poly[j] = nmod_sub(poly[j], poly[j + 1], mod);
    }
    else if (c != 0)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                NMOD_ADDMUL(poly[j], poly[j + 1], c, mod);
    }
}

void
nmod_poly_taylor_shift_horner(nmod_poly_t g, const nmod_poly_t f, ulong c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift_horner(g->coeffs, c, g->length, g->mod);
}
