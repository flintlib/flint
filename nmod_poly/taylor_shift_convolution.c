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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_taylor_shift_convolution(mp_ptr p, mp_limb_t c, len_t len, nmod_t mod)
{
    len_t i, n = len - 1;
    mp_limb_t f, d;
    mp_ptr t, u;

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
    mp_limb_t c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift_convolution(g->coeffs, c, g->length, g->mod);
}
