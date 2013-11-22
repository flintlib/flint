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
#include "long_extras.h"

void
_nmod_poly_powmod_x_ui_preinv (mp_ptr res, ulong e, mp_srcptr f, slong lenf,
                               mp_srcptr finv, slong lenfinv, nmod_t mod)
{
    mp_ptr T, Q;
    slong lenT, lenQ, window;
    int i, l, c;

    lenT = 2 * lenf - 3;
    lenQ = FLINT_MAX(lenT - lenf + 1, 1);

    T = _nmod_vec_init(lenT + lenQ);
    Q = T + lenT;

    flint_mpn_zero (res, lenf - 1);
    res[0] = WORD(1);

    l = (int) z_sizeinbase (lenf - 1, 2) - 2;
    window = WORD(0);
    window = (WORD(1) << l);
    c = l;
    i = (int) FLINT_BIT_COUNT(e) - 2;
    if (i <= l)
    {
      window = WORD(0);
      window = (WORD(1) << i);
      c = i;
      l = i;
    }

    if (c == 0)
    {
        _nmod_poly_shift_left(T, res, lenf - 1, window);
        _nmod_poly_divrem_newton_n_preinv(Q, res, T, lenf - 1 + window, f,
                                          lenf, finv, lenfinv, mod);
        c = l + 1;
        window= WORD(0);
    }

    for (; i >= 0; i--)
    {
        _nmod_poly_mul(T, res, lenf - 1, res, lenf - 1, mod);
        _nmod_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f,
                                          lenf, finv, lenfinv, mod);

        c--;
        if (e & (UWORD(1) << i))
        {
            if (window == WORD(0) && i <= l - 1)
                c = i;
            if ( c >= 0)
              window = window | (WORD(1) << c);
        }
        else if (window == WORD(0))
            c= l + 1;
        if (c == 0)
        {
            _nmod_poly_shift_left(T, res, lenf - 1, window);
            _nmod_poly_divrem_newton_n_preinv(Q, res, T, lenf - 1 + window, f,
                                              lenf, finv, lenfinv, mod);

            c= l + 1;
            window= WORD(0);
        }
    }

    _nmod_vec_clear(T);
}


void
nmod_poly_powmod_x_ui_preinv(nmod_poly_t res, ulong e, const nmod_poly_t f,
                             const nmod_poly_t finv)
{
    slong lenf = f->length;
    slong trunc = lenf - 1;
    nmod_poly_t tmp;

    if (lenf == 0)
    {
        flint_printf("Exception (nmod_poly_powmod_x_ui_preinv). Divide by zero.\n");
        abort();
    }

    if (lenf == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (lenf == 2)
    {
        nmod_poly_t r, poly;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_init_preinv(r, res->mod.n, res->mod.ninv);
        nmod_poly_init2_preinv(poly, res->mod.n, res->mod.ninv, 2);
        nmod_poly_set_coeff_ui (poly, 1, 1);
        nmod_poly_divrem(tmp, r, poly, f);
        nmod_poly_powmod_ui_binexp_preinv(res, r, e, f, finv);
        nmod_poly_clear(tmp);
        nmod_poly_clear(r);
        nmod_poly_clear(poly);
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
            nmod_poly_t r;
            nmod_poly_init2_preinv (r, res->mod.n, res->mod.ninv, 2);
            nmod_poly_set_coeff_ui (r, 1, 1);
            nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
            nmod_poly_divrem (tmp, res, r, f);
            nmod_poly_clear(tmp);
            nmod_poly_clear(r);
        }
        else
        {
            nmod_poly_init2_preinv (tmp, res->mod.n, res->mod.ninv, 3);
            nmod_poly_set_coeff_ui (tmp, 1, 1);
            nmod_poly_mulmod (res, tmp, tmp, f);
            nmod_poly_clear(tmp);
        }
        return;
    }

    if ((res == f) || (res == finv))
    {
        nmod_poly_init2(tmp, res->mod.n, trunc);
        _nmod_poly_powmod_x_ui_preinv(tmp->coeffs, e, f->coeffs, lenf,
                                      finv->coeffs, finv->length, f->mod);
        nmod_poly_swap(res, tmp);
        nmod_poly_clear(tmp);
    }
    else
    {
        nmod_poly_fit_length(res, trunc);
        _nmod_poly_powmod_x_ui_preinv(res->coeffs, e, f->coeffs, lenf,
                                      finv->coeffs, finv->length, f->mod);
    }

    res->length = trunc;
    _nmod_poly_normalise(res);
}
