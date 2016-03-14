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
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
void
_TEMPLATE(T, poly_compose_mod) (TEMPLATE(T, struct) * res,
                                const TEMPLATE(T, struct) * f, slong lenf,
                                const TEMPLATE(T, struct) * g,
                                const TEMPLATE(T, struct) * h, slong lenh,
                                const TEMPLATE(T, ctx_t) ctx)
{
    if (lenh < TEMPLATE(CAP_T, COMPOSE_MOD_LENH_CUTOFF) || lenf >= lenh)
        _TEMPLATE(T, poly_compose_mod_horner) (res, f, lenf, g, h, lenh, ctx);
    else
        _TEMPLATE(T, poly_compose_mod_brent_kung) (res, f, lenf, g, h, lenh,
                                                   ctx);
}

void
TEMPLATE(T, poly_compose_mod) (TEMPLATE(T, poly_t) res,
                               const TEMPLATE(T, poly_t) poly1,
                               const TEMPLATE(T, poly_t) poly2,
                               const TEMPLATE(T, poly_t) poly3,
                               const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) inv3;
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    TEMPLATE(T, struct) * ptr2;

    if (len3 == 0)
    {
        TEMPLATE_PRINTF("Exception: division by zero in %s_poly_compose_mod\n",
                        T);
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    if (len1 == 1)
    {
        TEMPLATE(T, poly_set) (res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        TEMPLATE(T, poly_t) tmp;
        TEMPLATE(T, poly_init) (tmp, ctx);
        TEMPLATE(T, poly_compose_mod) (tmp, poly1, poly2, poly3, ctx);
        TEMPLATE(T, poly_swap) (tmp, res, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
        return;
    }

    ptr2 = _TEMPLATE(T, vec_init) (vec_len, ctx);

    if (len2 <= len)
    {
        _TEMPLATE(T, vec_set) (ptr2, poly2->coeffs, len2, ctx);
        _TEMPLATE(T, vec_zero) (ptr2 + len2, len - len2, ctx);
    }
    else
    {
        TEMPLATE(T, init) (inv3, ctx);
        TEMPLATE(T, inv) (inv3, poly3->coeffs + len, ctx);
        _TEMPLATE(T, poly_rem) (ptr2, poly2->coeffs, len2,
                                poly3->coeffs, len3, inv3, ctx);
        TEMPLATE(T, clear) (inv3, ctx);
    }

    TEMPLATE(T, poly_fit_length) (res, len, ctx);
    _TEMPLATE(T, poly_compose_mod) (res->coeffs,
                                    poly1->coeffs, len1, ptr2, poly3->coeffs,
                                    len3, ctx);
    _TEMPLATE(T, poly_set_length) (res, len, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);

    _TEMPLATE(T, vec_clear) (ptr2, vec_len, ctx);
}


#endif
