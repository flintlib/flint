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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fq_poly.h"

void
_fq_poly_compose_mod(fq_struct * res, 
                     const fq_struct * f, slong lenf, 
                     const fq_struct * g,
                     const fq_struct * h, slong lenh, 
                     const fq_ctx_t ctx)
{
    if (lenh < FQ_COMPOSE_MOD_LENH_CUTOFF || lenf >= lenh)
        _fq_poly_compose_mod_horner(res, f, lenf, g, h, lenh, ctx);
    else
        _fq_poly_compose_mod_brent_kung(res, f, lenf, g, h, lenh, ctx);
}

void
fq_poly_compose_mod(fq_poly_t res, const fq_poly_t poly1,
                    const fq_poly_t poly2, const fq_poly_t poly3,
                    const fq_ctx_t ctx)
{
    fq_t inv3;
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    fq_struct * ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception: division by zero in fq_poly_compose_mod\n");
        abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        fq_poly_zero(res, ctx);
        return;
    }

    if (len1 == 1)
    {
        fq_poly_set(res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        fq_poly_t tmp;
        fq_poly_init(tmp, ctx);
        fq_poly_compose_mod(tmp, poly1, poly2, poly3, ctx);
        fq_poly_swap(tmp, res, ctx);
        fq_poly_clear(tmp, ctx);
        return;
    }

    ptr2 = _fq_vec_init(vec_len, ctx);

    if (len2 <= len)
    {
        _fq_vec_set(ptr2, poly2->coeffs, len2, ctx);
        _fq_vec_zero(ptr2 + len2, len - len2, ctx);
    }
    else
    {
        fq_init(inv3, ctx);
        fq_inv(inv3, poly3->coeffs + len, ctx);
        _fq_poly_rem(ptr2, poly2->coeffs, len2,
                               poly3->coeffs, len3, inv3, ctx);
        fq_clear(inv3, ctx);
    }

    fq_poly_fit_length(res, len, ctx);
    _fq_poly_compose_mod(res->coeffs,
                                   poly1->coeffs, len1, 
                                   ptr2, 
                                   poly3->coeffs, len3,
                                   ctx);
    _fq_poly_set_length(res, len, ctx);
    _fq_poly_normalise(res, ctx);

    _fq_vec_clear(ptr2, vec_len, ctx);
}
