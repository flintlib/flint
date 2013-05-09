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
   Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_compose_series(mp_ptr res, mp_srcptr poly1, len_t len1, 
                            mp_srcptr poly2, len_t len2, len_t n, nmod_t mod)
{
    if (len1 < 24 || len2 < 8)
        _nmod_poly_compose_series_horner(res, poly1, len1,
                                                    poly2, len2, n, mod);
    else
        _nmod_poly_compose_series_brent_kung(res, poly1, len1,
                                                    poly2, len2, n, mod);
}

void
nmod_poly_compose_series(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, len_t n)
{
    len_t len1 = poly1->length;
    len_t len2 = poly2->length;
    len_t lenr;

    if (len2 != 0 && poly2->coeffs[0] != 0)
    {
        printf("Exception (nmod_poly_compose_series). Inner polynomial \n"
               "must have zero constant term.\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = 1;
        _nmod_poly_normalise(res);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        nmod_poly_fit_length(res, lenr);
        _nmod_poly_compose_series(res->coeffs, poly1->coeffs, len1, 
                                        poly2->coeffs, len2, lenr, res->mod);
        res->length = lenr;
        _nmod_poly_normalise(res);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2_preinv(t, res->mod.n, res->mod.ninv, lenr);
        _nmod_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, res->mod);
        t->length = lenr;
        _nmod_poly_normalise(t);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
}
