/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

#define __nmod_poly_get_coeff(p,i) ((p)->coeffs[(i)])
#define __fmpz_mod_poly_get_coeff(p,i) ((p)->coeffs + (i))


void TEMPLATE(T, embed_mono_to_dual_matrix)(TEMPLATE(B, mat_t) res,
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

void TEMPLATE(T, embed_dual_to_mono_matrix)(TEMPLATE(B, mat_t) res,
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
    TEMPLATE(T, embed_mul_matrix)(mul_mat, d_ctx_inv, ctx);
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
void TEMPLATE(T, embed_trace_matrix)(TEMPLATE(B, mat_t) res,
                               const TEMPLATE(B, mat_t) basis,
                               const TEMPLATE(T, ctx_t) sub_ctx,
                               const TEMPLATE(T, ctx_t) sup_ctx)
{
    slong m = TEMPLATE(B, mat_ncols)(basis),
          n = TEMPLATE(B, mat_nrows)(basis);
    const TEMPLATE(B, poly_struct) *modulus = TEMPLATE(T, ctx_modulus)(sub_ctx);
    TEMPLATE(B, mat_t) m2d, d2m, tmp;

    TEMPLATE(B, mat_init)(m2d, n, n, TEMPLATE(B, poly_modulus)(modulus));
    TEMPLATE(B, mat_init)(d2m, m, m, TEMPLATE(B, poly_modulus)(modulus));
    TEMPLATE(B, mat_init)(tmp, m, n, TEMPLATE(B, poly_modulus)(modulus));

    TEMPLATE(T, embed_mono_to_dual_matrix)(m2d, sup_ctx);
    TEMPLATE(B, mat_transpose)(res, basis);
    TEMPLATE(T, embed_dual_to_mono_matrix)(d2m, sub_ctx);

    TEMPLATE(B, mat_mul)(tmp, res, m2d);
    TEMPLATE(B, mat_mul)(res, d2m, tmp);

    TEMPLATE(B, mat_clear)(m2d);
    TEMPLATE(B, mat_clear)(d2m);
    TEMPLATE(B, mat_clear)(tmp);
}

static inline
int __fmpz_mod_inv_degree(fmpz_t invd, slong d, const fmpz_t p)
{
    fmpz_set_si(invd, d);
    return fmpz_invmod(invd, invd, p);
}

static inline
int __nmod_inv_degree(fmpz_t invd, slong d, mp_limb_t p)
{
    mp_limb_t ud = d % p;
    if (!ud)
        return 0;
    ud = n_invmod(ud, p);
    fmpz_set_ui(invd, ud);
    return 1;
}

#define fmpz_mod_mat_entry_is_zero(mat, i, j) (fmpz_is_zero(fmpz_mod_mat_entry((mat), (i), (j))))
#define nmod_mat_entry_is_zero(mat, i, j) (nmod_mat_entry((mat), (i), (j)) == 0)

void TEMPLATE(T, embed_matrices)(TEMPLATE(B, mat_t) embed,
                                 TEMPLATE(B, mat_t) project,
                                 const TEMPLATE(T, t) gen_sub,
                                 const TEMPLATE(T, ctx_t) sub_ctx,
                                 const TEMPLATE(T, t) gen_sup,
                                 const TEMPLATE(T, ctx_t) sup_ctx,
                                 const TEMPLATE(B, poly_t) gen_minpoly)
{
    slong m = TEMPLATE(T, ctx_degree)(sub_ctx),
        n = TEMPLATE(T, ctx_degree)(sup_ctx),
        d = n / m;
    fmpz_t invd;
    TEMPLATE(T, ctx_t) gen_ctx;
    TEMPLATE(B, poly_t) gen_minpoly_cp;
    TEMPLATE(B, mat_t) gen2sub, sub2gen, gen2sup, sup2gen;

    /* Is there any good reason why the modulus argument to
       fq_ctx_init_modulus is not const? */
    TEMPLATE(B, poly_init)(gen_minpoly_cp, TEMPLATE(B, poly_modulus)(gen_minpoly));
    TEMPLATE(B, poly_set)(gen_minpoly_cp, gen_minpoly);
    fmpz_init(invd);
    TEMPLATE(T, ctx_init_modulus)(gen_ctx, gen_minpoly_cp, "x");
    TEMPLATE(B, mat_init)(gen2sub, m, m, TEMPLATE(B, poly_modulus)(gen_minpoly));
    TEMPLATE(B, mat_init)(sub2gen, m, m, TEMPLATE(B, poly_modulus)(gen_minpoly));
    TEMPLATE(B, mat_init)(gen2sup, n, m, TEMPLATE(B, poly_modulus)(gen_minpoly));
    TEMPLATE(B, mat_init)(sup2gen, m, n, TEMPLATE(B, poly_modulus)(gen_minpoly));

    /* Gen -> Sub */
    TEMPLATE(T, embed_composition_matrix)(gen2sub, gen_sub, sub_ctx);
    /* Sub -> Gen */
    TEMPLATE(T, embed_trace_matrix)(sub2gen, gen2sub, gen_ctx, sub_ctx);
    /* Gen -> Sup */
    TEMPLATE(T, embed_composition_matrix_sub)(gen2sup, gen_sup, sup_ctx, m);
    /* Sup -> Gen (trace) */
    TEMPLATE(T, embed_trace_matrix)(sup2gen, gen2sup, gen_ctx, sup_ctx);

    /* Correct the projection Sup -> Gen */
    /* If this is an isomorphism, there is no need for correction */
    if (d == 1) {}
    /* If the extension degree is invertible mod p, multiply trace by 1/d */
    else if (__TEMPLATE(B, inv_degree)(invd, d, TEMPLATE(B, poly_modulus)(gen_minpoly)))
    {
        TEMPLATE(B, mat_scalar_mul_fmpz)(sup2gen, sup2gen, invd);
    }
    /* Otherwise pre-multiply by an element of trace equal to 1 */
    else
    {
        int i;
        TEMPLATE(T, t) mul, trace;
        TEMPLATE(B, mat_t) column, tvec, mat_mul, tmp;

        TEMPLATE(T, init)(mul, sup_ctx);
        TEMPLATE(T, init)(trace, sup_ctx);
        TEMPLATE(B, mat_init)(tvec, n, 1, TEMPLATE(B, poly_modulus)(gen_minpoly));
        TEMPLATE(B, mat_init)(mat_mul, n, n, TEMPLATE(B, poly_modulus)(gen_minpoly));
        TEMPLATE(B, mat_init)(tmp, m, n, TEMPLATE(B, poly_modulus)(gen_minpoly));

        /* Look for a non-zero column in sup2gen
           (we know it has full rank, so first row is non-null) */
        for (i = 1; i < n; i++)
        {
            if (!TEMPLATE(B, mat_entry_is_zero(sup2gen, 0, i)))
                break;
        }

        /* Set mul to x^i */
        TEMPLATE(T, gen)(mul, sup_ctx);
        TEMPLATE(T, pow_ui)(mul, mul, i, sup_ctx);
        /* Set trace to its trace */
        TEMPLATE(B, mat_window_init)(column, sup2gen, 0, i, m, i+1);
        TEMPLATE(B, mat_mul)(tvec, gen2sup, column);
        TEMPLATE4(T, set, B, mat)(trace, tvec, sup_ctx);
        /* Get an element of trace 1 */
        TEMPLATE(T, div)(mul, mul, trace, sup_ctx);

        /* Correct the matrix */
        TEMPLATE(T, embed_mul_matrix)(mat_mul, mul, sup_ctx);
        TEMPLATE(B, mat_mul)(tmp, sup2gen, mat_mul);
        TEMPLATE(B, mat_swap)(tmp, sup2gen);

        TEMPLATE(B, mat_clear)(tmp);
        TEMPLATE(B, mat_clear)(mat_mul);
        TEMPLATE(B, mat_clear)(tvec);
        TEMPLATE(B, mat_window_clear)(column);
        TEMPLATE(T, clear)(mul, sup_ctx);
        TEMPLATE(T, clear)(trace, sup_ctx);
    }

    TEMPLATE(B, mat_mul)(embed, gen2sup, sub2gen);
    TEMPLATE(B, mat_mul)(project, gen2sub, sup2gen);

    TEMPLATE(B, mat_clear)(gen2sub);
    TEMPLATE(B, mat_clear)(sub2gen);
    TEMPLATE(B, mat_clear)(gen2sup);
    TEMPLATE(B, mat_clear)(sup2gen);
    TEMPLATE(T, ctx_clear)(gen_ctx);
    fmpz_clear(invd);
    TEMPLATE(B, poly_clear)(gen_minpoly_cp);
}


#endif
#endif
