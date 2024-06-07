/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"
#include "acb.h"
#include "acb_poly.h"

/* trivial test code: stress tests TBD */
TEST_FUNCTION_START(gr_poly_refine_roots, state)
{
    gr_poly_t f, g;
    acb_ptr r, z, z2, w;
    gr_ctx_t ctx;
    slong i, which;
    mag_t m;

    gr_ctx_init_complex_acb(ctx, 256);
    gr_poly_init(f, ctx);
    gr_poly_init(g, ctx);

    r = _acb_vec_init(3);
    z = _acb_vec_init(3);
    z2 = _acb_vec_init(3);
    w = _acb_vec_init(3);
    mag_init(m);

    for (i = 0; i < 4; i++)
        GR_IGNORE(gr_poly_set_coeff_si(f, i, i + 1, ctx));
    GR_IGNORE(gr_poly_derivative(g, f, ctx));

    for (which = 0; which < 4; which++)
    {
        acb_set_d_d(r + 0, -0.60582958618826802099, 0.0);
        acb_set_d_d(r + 1, -0.0720852069058659895, 0.63832673514837645798);
        acb_set_d_d(r + 2, -0.0720852069058659895, -0.63832673514837645798);

        acb_set_d_d(z + 0, -0.60582958618826802099 + 1e-4, 0.0);
        acb_set_d_d(z + 1, -0.0720852069058659895 + 1e-4, 0.63832673514837645798 + 1e-4);
        acb_set_d_d(z + 2, -0.0720852069058659895 - 2e-4, -0.63832673514837645798 + 0.5e-4);

        if (which == 0)
            GR_MUST_SUCCEED(_gr_poly_refine_roots_wdk(w, f->coeffs, 3, z, 0, ctx));
        else if (which == 1)
            GR_MUST_SUCCEED(_gr_poly_refine_roots_wdk(w, f->coeffs, 3, z, 1, ctx));
        else if (which == 2)
            GR_MUST_SUCCEED(_gr_poly_refine_roots_aberth(w, f->coeffs, g->coeffs, 3, z, 0, ctx));
        else
            GR_MUST_SUCCEED(_gr_poly_refine_roots_aberth(w, f->coeffs, g->coeffs, 3, z, 1, ctx));

        GR_MUST_SUCCEED(_gr_vec_sub(z2, z, w, 3, ctx));

        for (i = 0; i < 3; i++)
        {
            acb_sub(z + i, z2 + i, r + i, 53);
            acb_get_mag(m, z + i);

            if (mag_cmp_2exp_si(m, -20) > 0)
            {
                flint_printf("FAIL:\n");
                _gr_vec_print(r, 3, ctx); flint_printf("\n\n");
                _gr_vec_print(z2, 3, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }
    }

    _acb_vec_clear(r, 3);
    _acb_vec_clear(z, 3);
    _acb_vec_clear(z2, 3);
    _acb_vec_clear(w, 3);
    mag_clear(m);

    gr_poly_clear(f, ctx);
    gr_poly_clear(g, ctx);

    gr_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
