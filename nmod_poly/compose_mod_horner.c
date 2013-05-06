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

   Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_compose_mod_horner(mp_ptr res,
    mp_srcptr f, long lenf, mp_srcptr g, mp_srcptr h, long lenh, nmod_t mod)
{
    long i, len;
    mp_ptr t;

    if (lenh == 1)
        return;

    if (lenf == 1)
    {
        res[0] = f[0];
        return;
    }

    if (lenh == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(f, lenf, g[0], mod);
        return;
    }

    len = lenh - 1;
    i = lenf - 1;
    t = _nmod_vec_init(len);

    _nmod_vec_scalar_mul_nmod(res, g, len, f[i], mod);
    i--;
    if (i >= 0)
        res[0] = nmod_add(res[0], f[i], mod);

    while (i > 0)
    {
        i--;
        _nmod_poly_mulmod(t, res, len, g, len, h, lenh, mod);
        _nmod_poly_add(res, t, len, f + i, 1, mod);
    }

    _nmod_vec_clear(t);
}

void
nmod_poly_compose_mod_horner(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long len3 = poly3->length;
    long len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        printf("Exception (nmod_poly_compose_mod_horner). Division by zero.\n");
        abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod_horner(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_horner(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
