/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

#define __nmod_poly_get_coeff(p,i) ((p)->coeffs[(i)])
#define __fmpz_mod_poly_get_coeff(p,i) ((p)->coeffs + (i))


void TEMPLATE(T, mono_to_dual_matrix)(TEMPLATE(B, mat_t) res,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, n = TEMPLATE(T, ctx_degree)(ctx);
    TEMPLATE(B, poly_t) ctx_inv_rev, d_ctx;
    const TEMPLATE(B, poly_struct) *modulus = TEMPLATE(T, ctx_modulus)(ctx);

    TEMPLATE(B, poly_init)(ctx_inv_rev, TEMPLATE(B, poly_modulus)(modulus));
    TEMPLATE(B, poly_init)(d_ctx, TEMPLATE(B, poly_modulus)(modulus));

    /* Half of this is precomputed in ctx. Maybe a call to some
       internal Newton stuff could be enough to double it. */
    TEMPLATE(T, modulus_pow_series_inv)(ctx_inv_rev, ctx, 2*n);
    TEMPLATE(B, poly_derivative)(d_ctx, modulus);
    TEMPLATE(B, poly_reverse)(d_ctx, d_ctx, n);
    TEMPLATE(B, poly_mullow)(ctx_inv_rev, ctx_inv_rev, d_ctx, 2*n);

    TEMPLATE(B, mat_zero)(res);
    for (i = 0; i < n; i++)
        for (j = 0; j < n && i+j < ctx_inv_rev->length; j++)
            TEMPLATE(B, mat_set_entry)(res, i, j,
                                       __TEMPLATE(B, poly_get_coeff)(ctx_inv_rev, i + j));

    TEMPLATE(B, poly_clear)(ctx_inv_rev);
    TEMPLATE(B, poly_clear)(d_ctx);
}

void TEMPLATE(T, dual_to_mono_matrix)(TEMPLATE(B, mat_t) res,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, n = TEMPLATE(T, ctx_degree)(ctx);
    TEMPLATE(T, t) d_ctx, d_ctx_inv;
    const TEMPLATE(B, poly_struct) *modulus = TEMPLATE(T, ctx_modulus)(ctx);
    TEMPLATE(B, mat_t) mul_mat, tmp;

    TEMPLATE(T, init)(d_ctx, ctx);
    TEMPLATE(T, init)(d_ctx_inv, ctx);
    TEMPLATE(B, mat_init)(mul_mat, n, n, TEMPLATE(B, poly_modulus)(modulus));
    TEMPLATE(B, mat_init)(tmp, n, n, TEMPLATE(B, poly_modulus)(modulus));

    TEMPLATE(B, mat_zero)(tmp);
    for (i = 0; i < n; i++)
        for (j = 0; j < n - i; j++)
            TEMPLATE(B, mat_set_entry)(tmp, i, j,
                                       __TEMPLATE(B, poly_get_coeff)(modulus, i + j + 1));

    TEMPLATE(T, modulus_derivative_inv)(d_ctx, d_ctx_inv, ctx);
    TEMPLATE(T, mul_matrix)(mul_mat, d_ctx_inv, ctx);
    TEMPLATE(B, mat_mul)(res, mul_mat, tmp);

    TEMPLATE(T, clear)(d_ctx, ctx);
    TEMPLATE(T, clear)(d_ctx_inv, ctx);
    TEMPLATE(B, mat_clear)(mul_mat);
    TEMPLATE(B, mat_clear)(tmp);
}

/* Given:

   - a context `sup_ctx` defining a field K of degree n,
   - a context `sub_ctx` defining a subfield k of degree m|n,
   - the n×m change-of-basis matrix from k to K,

   Return the m×n matrix of the trace from K to k.
   If m=n, res is the inverse of basis. */
void TEMPLATE(T, trace_matrix)(TEMPLATE(B, mat_t) res,
                               const TEMPLATE(B, mat_t) basis,
                               const TEMPLATE(T, ctx_t) sub_ctx,
                               const TEMPLATE(T, ctx_t) sup_ctx)
{
    slong m = TEMPLATE(B, mat_ncols)(basis), 
        n = TEMPLATE(B, mat_nrows)(basis);
    const TEMPLATE(B, poly_struct) *modulus = TEMPLATE(T, ctx_modulus)(sub_ctx);
    TEMPLATE(B, mat_t) m2d, d2m;

    TEMPLATE(B, mat_init)(m2d, n, n, TEMPLATE(B, poly_modulus)(modulus));
    TEMPLATE(B, mat_init)(d2m, m, m, TEMPLATE(B, poly_modulus)(modulus));

    TEMPLATE(T, mono_to_dual_matrix)(m2d, sup_ctx);
    TEMPLATE(B, mat_transpose)(res, basis);
    TEMPLATE(T, dual_to_mono_matrix)(d2m, sub_ctx);

    TEMPLATE(B, mat_mul)(res, res, m2d);
    TEMPLATE(B, mat_mul)(res, d2m, res);

    TEMPLATE(B, mat_clear)(m2d);
    TEMPLATE(B, mat_clear)(d2m);
}


#endif
#endif
