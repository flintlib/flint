/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "nmod_poly_factor_gr.h"

/*
    Helper function for finding roots. The roots of a monic f are written
    with exponent given in mult to r. Uses Rabin's Las Vegas algorithm
    via gcd computations with (x + delta)^((p-1)/2) - 1.
*/
void nmod_poly_roots(nmod_poly_factor_t r, const nmod_poly_t f,
                                                         int with_multiplicity)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_vec_t roots;
    fmpz_vec_t mult;
    nmod_t mod = f->mod;
    slong i, num;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, f);
    gr_vec_init(roots, 0, ctx);
    fmpz_vec_init(mult, 0);

    GR_MUST_SUCCEED(gr_poly_roots_finite_field(roots, mult, P,
        0, ctx));

    num = roots->length;
    r->num = 0;
    nmod_poly_factor_fit_length(r, num);

    for (i = 0; i < num; i++)
    {
        /* the linear factor x - root */
        nmod_poly_struct * p = r->p + i;
        nmod_poly_fit_length(p, 2);
        p->coeffs[0] = nmod_neg(((nn_srcptr) roots->entries)[i], mod);
        p->coeffs[1] = 1;
        p->length = 2;
        r->exp[i] = with_multiplicity ? fmpz_get_si(mult->entries + i) : 1;
    }

    r->num = num;

    gr_vec_clear(roots, ctx);
    fmpz_vec_clear(mult);
    gr_ctx_clear(ctx);
}
