/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_embed.h"
#include "fq_poly.h"


void fq_embed_mono_to_dual_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)
{
    slong i, j, n = fq_ctx_degree(ctx);
    fmpz_mod_poly_t ctx_inv_rev, d_ctx;
    const fmpz_mod_poly_struct *modulus = fq_ctx_modulus(ctx);

    fmpz_mod_poly_init(ctx_inv_rev, ctx->ctxp);
    fmpz_mod_poly_init(d_ctx, ctx->ctxp);

    /* Half of this is precomputed in ctx. Maybe a call to some
       internal Newton stuff could be enough to double it. */
    fq_modulus_pow_series_inv(ctx_inv_rev, ctx, 2*n);
    fmpz_mod_poly_derivative(d_ctx, modulus, ctx->ctxp);
    fmpz_mod_poly_reverse(d_ctx, d_ctx, n, ctx->ctxp);
    fmpz_mod_poly_mullow(ctx_inv_rev, ctx_inv_rev, d_ctx, 2*n, ctx->ctxp);

    fmpz_mod_mat_zero(res);
    for (i = 0; i < n; i++)
        for (j = 0; j < n && i+j < ctx_inv_rev->length; j++)
            fmpz_mod_mat_set_entry(res, i, j, ctx_inv_rev->coeffs + i + j);

    fmpz_mod_poly_clear(ctx_inv_rev, ctx->ctxp);
    fmpz_mod_poly_clear(d_ctx, ctx->ctxp);
}

void fq_embed_dual_to_mono_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)
{
    slong i, j, n = fq_ctx_degree(ctx);
    fq_t d_ctx, d_ctx_inv;
    const fmpz_mod_poly_struct *modulus = fq_ctx_modulus(ctx);
    fmpz_mod_mat_t mul_mat, tmp;

    fq_init(d_ctx, ctx);
    fq_init(d_ctx_inv, ctx);
    fmpz_mod_mat_init(mul_mat, n, n, fq_ctx_prime(ctx));
    fmpz_mod_mat_init(tmp, n, n, fq_ctx_prime(ctx));

    fmpz_mod_mat_zero(tmp);
    for (i = 0; i < n; i++)
        for (j = 0; j < n - i; j++)
            fmpz_mod_mat_set_entry(tmp, i, j, modulus->coeffs + i + j + 1);

    fq_modulus_derivative_inv(d_ctx, d_ctx_inv, ctx);
    fq_embed_mul_matrix(mul_mat, d_ctx_inv, ctx);
    fmpz_mod_mat_mul(res, mul_mat, tmp);

    fq_clear(d_ctx, ctx);
    fq_clear(d_ctx_inv, ctx);
    fmpz_mod_mat_clear(mul_mat);
    fmpz_mod_mat_clear(tmp);
}

/* Given:

   - a context `sup_ctx` defining a field K of degree n,
   - a context `sub_ctx` defining a subfield k of degree m|n,
   - the n×m change-of-basis matrix from k to K,

   Return the m×n matrix of the trace from K to k.
   If m=n, res is the inverse of basis. */
void fq_embed_trace_matrix(fmpz_mod_mat_t res,
                               const fmpz_mod_mat_t basis,
                               const fq_ctx_t sub_ctx,
                               const fq_ctx_t sup_ctx)
{
    slong m = fmpz_mod_mat_ncols(basis); 
    slong n = fmpz_mod_mat_nrows(basis);
    fmpz_mod_mat_t m2d, d2m, tmp;

    fmpz_mod_mat_init(m2d, n, n, fmpz_mod_ctx_modulus(sub_ctx->ctxp));
    fmpz_mod_mat_init(d2m, m, m, fmpz_mod_ctx_modulus(sub_ctx->ctxp));
    fmpz_mod_mat_init(tmp, m, n, fmpz_mod_ctx_modulus(sub_ctx->ctxp));

    fq_embed_mono_to_dual_matrix(m2d, sup_ctx);
    fmpz_mod_mat_transpose(res, basis);
    fq_embed_dual_to_mono_matrix(d2m, sub_ctx);

    fmpz_mod_mat_mul(tmp, res, m2d);
    fmpz_mod_mat_mul(res, d2m, tmp);

    fmpz_mod_mat_clear(m2d);
    fmpz_mod_mat_clear(d2m);
    fmpz_mod_mat_clear(tmp);
}

void fq_embed_matrices(fmpz_mod_mat_t embed,
                                 fmpz_mod_mat_t project,
                                 const fq_t gen_sub,
                                 const fq_ctx_t sub_ctx,
                                 const fq_t gen_sup,
                                 const fq_ctx_t sup_ctx,
                                 const fmpz_mod_poly_t gen_minpoly)
{
    const fmpz_mod_ctx_struct * ctxp = sub_ctx->ctxp;
    slong m = fq_ctx_degree(sub_ctx);
    slong n = fq_ctx_degree(sup_ctx);
    slong d = n / m;
    fmpz_t invd;
    fq_ctx_t gen_ctx;
    fmpz_mod_poly_t gen_minpoly_cp;
    fmpz_mod_mat_t gen2sub, sub2gen, gen2sup, sup2gen;

    /* Is there any good reason why the modulus argument to
       fq_ctx_init_modulus is not const? */
    fmpz_mod_poly_init(gen_minpoly_cp, ctxp);
    fmpz_mod_poly_set(gen_minpoly_cp, gen_minpoly, ctxp);
    fmpz_init(invd);
    fq_ctx_init_modulus(gen_ctx, gen_minpoly_cp, ctxp, "x");
    fmpz_mod_mat_init(gen2sub, m, m, fmpz_mod_ctx_modulus(ctxp));
    fmpz_mod_mat_init(sub2gen, m, m, fmpz_mod_ctx_modulus(ctxp));
    fmpz_mod_mat_init(gen2sup, n, m, fmpz_mod_ctx_modulus(ctxp));
    fmpz_mod_mat_init(sup2gen, m, n, fmpz_mod_ctx_modulus(ctxp));

    /* Gen -> Sub */
    fq_embed_composition_matrix(gen2sub, gen_sub, sub_ctx);
    /* Sub -> Gen */
    fq_embed_trace_matrix(sub2gen, gen2sub, gen_ctx, sub_ctx);
    /* Gen -> Sup */
    fq_embed_composition_matrix_sub(gen2sup, gen_sup, sup_ctx, m);
    /* Sup -> Gen (trace) */
    fq_embed_trace_matrix(sup2gen, gen2sup, gen_ctx, sup_ctx);

    /* Correct the projection Sup -> Gen */
    /* If this is an isomorphism, there is no need for correction */
    if (d == 1) {}
    /* If the extension degree is invertible mod p, multiply trace by 1/d */
    else if (fmpz_set_si(invd, d),
             fmpz_invmod(invd, invd, fmpz_mod_ctx_modulus(ctxp)))
    {
        fmpz_mod_mat_scalar_mul_fmpz(sup2gen, sup2gen, invd);
    }
    /* Otherwise pre-multiply by an element of trace equal to 1 */
    else
    {
        int i;
        fq_t mul, trace;
        fmpz_mod_mat_t column, tvec, mat_mul, tmp;
        
        fq_init(mul, sup_ctx);
        fq_init(trace, sup_ctx);
        fmpz_mod_mat_init(tvec, n, 1, fmpz_mod_ctx_modulus(ctxp));
        fmpz_mod_mat_init(mat_mul, n, n, fmpz_mod_ctx_modulus(ctxp));
        fmpz_mod_mat_init(tmp, m, n, fmpz_mod_ctx_modulus(ctxp));

        /* Look for a non-zero column in sup2gen
           (we know it has full rank, so first row is non-null) */
        for (i = 1; i < n; i++)
        {
            if (!fmpz_is_zero(fmpz_mod_mat_entry(sup2gen, 0, i)))
                break;
        }
        
        /* Set mul to x^i */
        fq_gen(mul, sup_ctx);
        fq_pow_ui(mul, mul, i, sup_ctx);
        /* Set trace to its trace */
        fmpz_mod_mat_window_init(column, sup2gen, 0, i, m, i+1);
        fmpz_mod_mat_mul(tvec, gen2sup, column);
        fq_set_fmpz_mod_mat(trace, tvec, sup_ctx);
        /* Get an element of trace 1 */
        fq_div(mul, mul, trace, sup_ctx);
        
        /* Correct the matrix */
        fq_embed_mul_matrix(mat_mul, mul, sup_ctx);
        fmpz_mod_mat_mul(tmp, sup2gen, mat_mul);
        fmpz_mod_mat_swap(tmp, sup2gen);
        
        fmpz_mod_mat_clear(tmp);
        fmpz_mod_mat_clear(mat_mul);
        fmpz_mod_mat_clear(tvec);
        fmpz_mod_mat_window_clear(column);
        fq_clear(mul, sup_ctx);
        fq_clear(trace, sup_ctx);
    }

    fmpz_mod_mat_mul(embed, gen2sup, sub2gen);
    fmpz_mod_mat_mul(project, gen2sub, sup2gen);

    fmpz_mod_mat_clear(gen2sub);
    fmpz_mod_mat_clear(sub2gen);
    fmpz_mod_mat_clear(gen2sup);
    fmpz_mod_mat_clear(sup2gen);
    fq_ctx_clear(gen_ctx);
    fmpz_clear(invd);
    fmpz_mod_poly_clear(gen_minpoly_cp, ctxp);
}

