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


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_iterated_frobenius_preinv) (TEMPLATE(T, poly_t) * rop,
                                             slong n,
                                             const TEMPLATE(T, poly_t) v,
                                             const TEMPLATE(T, poly_t) vinv,
                                             const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    fmpz_t q;
    TEMPLATE(T, mat_t) HH;

    fmpz_init(q);

    TEMPLATE(T, ctx_order) (q, ctx);
    TEMPLATE(T, poly_gen) (rop[0], ctx);

    if (TEMPLATE(CAP_T, POLY_ITERATED_FROBENIUS_CUTOFF) (ctx, v->length))
    {
        TEMPLATE(T, mat_init) (HH, n_sqrt(v->length - 1) + 1, v->length - 1,
                               ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (rop[1], rop[0], q, 0, v,
                                                      vinv, ctx);
        TEMPLATE(T, poly_precompute_matrix) (HH, rop[1], v, vinv, ctx);
        for (i = 2; i < n; i++)
        {
            TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)
                (rop[i], rop[i - 1], HH, v, vinv, ctx);

        }
        TEMPLATE(T, mat_clear) (HH, ctx);
    }
    else
    {
        for (i = 1; i < n; i++)
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (rop[i], rop[i - 1],
                                                          q, 0, v, vinv, ctx);

    }

    fmpz_clear(q);

}

#endif
