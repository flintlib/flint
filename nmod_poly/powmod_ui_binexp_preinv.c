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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

******************************************************************************/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_powmod_ui_binexp_preinv (mp_ptr res, mp_srcptr poly,
                                    ulong e, mp_srcptr f, slong lenf,
                                    mp_srcptr finv, slong lenfinv, nmod_t mod)
{
    mp_ptr T, Q;
    slong lenT, lenQ;
    int i;

    if (lenf == 2)
    {
        res[0] = n_powmod2_ui_preinv(poly[0], e, mod.n, mod.ninv);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = FLINT_MAX(lenT - lenf + 1, 1);

    T = _nmod_vec_init(lenT + lenQ);
    Q = T + lenT;

    _nmod_vec_set(res, poly, lenf - 1);

    for (i = ((int) FLINT_BIT_COUNT(e) - 2); i >= 0; i--)
    {
        _nmod_poly_mul(T, res, lenf - 1, res, lenf - 1, mod);
        _nmod_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f,
                                          lenf, finv, lenfinv, mod);

        if (e & (UWORD(1) << i))
        {
            _nmod_poly_mul(T, res, lenf - 1, poly, lenf - 1, mod);
            _nmod_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f,
                                              lenf, finv, lenfinv, mod);
        }
    }

    _nmod_vec_clear(T);
}


void
nmod_poly_powmod_ui_binexp_preinv(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e,
                           const nmod_poly_t f, const nmod_poly_t finv)
{
    mp_ptr p;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int pcopy = 0;

    if (lenf == 0)
    {
        flint_printf("Exception (nmod_poly_powmod_ui_binexp_preinv). Divide by zero.\n");
        abort();
    }

    if (len >= lenf)
    {
        nmod_poly_t t, r;
        nmod_poly_init_preinv(t, res->mod.n, res->mod.ninv);
        nmod_poly_init_preinv(r, res->mod.n, res->mod.ninv);
        nmod_poly_divrem(t, r, poly, f);
        nmod_poly_powmod_ui_binexp_preinv(res, r, e, f, finv);
        nmod_poly_clear(t);
        nmod_poly_clear(r);
        return;
    }

    if (e <= 2)
    {
        if (e == UWORD(0))
        {
            nmod_poly_fit_length(res, 1);
            res->coeffs[0] = UWORD(1);
            res->length = 1;
        }
        else if (e == UWORD(1))
        {
            nmod_poly_set(res, poly);
        }
        else
            nmod_poly_mulmod_preinv(res, poly, poly, f, finv);
        return;
    }

    if (lenf == 1 || len == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len < trunc)
    {
        p = _nmod_vec_init(trunc);
        flint_mpn_copyi(p, poly->coeffs, len);
        flint_mpn_zero(p + len, trunc - len);
        pcopy = 1;
    } else
        p = poly->coeffs;

    if ((res == poly && !pcopy) || (res == f) || (res == finv))
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly->mod.n, trunc);
        _nmod_poly_powmod_ui_binexp_preinv(t->coeffs,
            p, e, f->coeffs, lenf, finv->coeffs, finv->length, poly->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
    else
    {
        nmod_poly_fit_length(res, trunc);
        _nmod_poly_powmod_ui_binexp_preinv(res->coeffs,
            p, e, f->coeffs, lenf, finv->coeffs, finv->length, poly->mod);
    }

    if (pcopy)
        _nmod_vec_clear(p);

    res->length = trunc;
    _nmod_poly_normalise(res);
}
