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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_pow(mp_ptr res, mp_srcptr poly, long len, ulong e, nmod_t mod)
{
    _nmod_poly_pow_binexp(res, poly, len, e, mod);
}

void
nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e)
{
    const long len = poly->length;
    long rlen;

    if ((len < 2) | (e < 3UL))
    {
        if (len == 0)
            nmod_poly_zero(res);
        else if (len == 1)
        {
            nmod_poly_fit_length(res, 1);
            res->coeffs[0] = n_powmod2_ui_preinv(poly->coeffs[0], e,
                poly->mod.n, poly->mod.ninv);
            res->length = 1;
            _nmod_poly_normalise(res);
        }
        else if (e == 0UL)
        {
            nmod_poly_set_coeff_ui(res, 0, 1UL);
            res->length = 1;
            _nmod_poly_normalise(res);
        }
        else if (e == 1UL)
            nmod_poly_set(res, poly);
        else  /* e == 2UL */
            nmod_poly_mul(res, poly, poly);

        return;
    }

    rlen = (long) e * (len - 1) + 1;

    if (res != poly)
    {
        nmod_poly_fit_length(res, rlen);
        _nmod_poly_pow(res->coeffs, poly->coeffs, len, e, poly->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly->mod.n, rlen);
        _nmod_poly_pow(t->coeffs, poly->coeffs, len, e, poly->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }

    res->length = rlen;
    _nmod_poly_normalise(res);
}
