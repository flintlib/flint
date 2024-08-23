/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* #include "arb.h" */
/* #include "acb.h" */
/* #include "arb_mat.h" */
/* #include "acb_theta.h" */

/* #define ACB_THETA_QL_TRY 10 */

/* /\* Assume: */
/*    - vectors zs always starts with zero. *\/ */

/* int acb_theta_ql_setup(acb_ptr rts, acb_ptr ts, slong * t_indices, acb_srcptr zs, slong nb, */
/*     const acb_mat_t tau, slong nb_steps, int all, int sqr, slong prec) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     slong n = 1 << g; */
/*     acb_theta_ctx_t ctx; */
/*     acb_theta_ctx_t new_ctx; */
/*     slong guard; */
/*     int * z_done; */
/*     slong nb_not_done; */
/*     slong * not_done; */
/*     int all_done = 0; */
/*     slong j, k; */

/*     FLINT_ASSERT(_acb_vec_is_zero(zs, g)); */

/*     acb_theta_ctx_init(ctx, nb, g); */
/*     z_done = flint_malloc(nb * sizeof(int)); */
/*     not_done = flint_malloc(nb * sizeof(slong)); */

/*     for (guard = 16; guard < prec; guard *= 2) */
/*     { */
/*         /\* Set context at precision guard *\/ */
/*         acb_theta_ctx_set_tau(ctx, tau, guard); */
/*         for (j = 0; j < nb; j++) */
/*         { */
/*             acb_theta_ctx_set_z(ctx, zs + j * g, j, guard); */
/*             acb_theta_ctx_use_t(ctx, j) = 0; */
/*         } */

/*         /\* See whether t = 0 works for every z; if so, exit loop *\/ */
/*         acb_theta_ql_roots(rts, ctx, guard); */
/*         all_done = 1; */
/*         nb_not_done = 0; */
/*         for (j = 0; j < nb; j++) */
/*         { */
/*             z_done[j] = 1; */
/*             for (k = 0; (k < n * nb_steps) && z_done[j]; k++) */
/*             { */
/*                 z_done[j] = z_done[j] && !acb_contains_zero(&rts[j * n * nb_steps + k]); */
/*             } */
/*             if (!z_done[j]) */
/*             { */
/*                 all_done = 0; */
/*                 acb_theta_ctx_use_t(ctx, j) = 1; */
/*                 not_done[nb_not_done] = j; */
/*                 nb_not_done += 1; */
/*             } */
/*         } */
/*         if (all_done) */
/*         { */
/*             break; */
/*         } */

/*         /\* Otherwise, pick auxiliary vector t. *\/ */
/*         /\* We must include zero in the "not done" list. *\/ */
/*         if (z_done[0]) */
/*         { */
/*             not_done[nb_not_done] = 0; */
/*             nb_not_done += 1; */
/*         } */
/*         for (k = 0; (k < ACB_THETA_QL_TRY) && !all_done; k++) */
/*         { */
/*             for (j = 0; j < g; j++) */
/*             { */
/*                 arb_urandom(acb_realref(&t[j]), state, prec); */
/*             } */
/*             acb_theta_ctx_init(new_ctx, 3 * nb_not_done, g); */
/*             acb_theta_ctx_copy_tau(new_ctx, ctx); */
/*             for (j = 0; j < nb_not_done; j++) */
/*             { */
/*                 acb_theta_ctx_copy_z(new_ctx, 3 * j, not_done[j]); */
/*             } */
/*             acb_theta_ctx_set_t(new_ctx, t, guard); */
/*             acb_theta_ql_roots(rts, new_ctx, guard); */
/*             all_done = 1; */
/*             for (k = 0; (k < 2 * n * nb_steps * nb_not_done) && all_done; k++) */
/*             { */
/*                 all_done = all_done && !acb_contains_zero(&rts[k]); */
/*             } */
/*             /\* Copy roots *\/ */
/*         } */
/*     } */

/*     return all_done; */
/* } */

/* void */
/* acb_theta_ctx_set_z_ql(acb_theta_ctx_t ctx, acb_srcptr z, slong j, slong prec) */
/* { */
/*     slong g = acb_theta_ctx_g(ctx); */
/*     slong n = 1 << g; */
/*     slong lp = ACB_THETA_LOW_PREC; */
/*     int z_is_zero = _acb_vec_is_zero(z, g); */
/*     int z_is_real = _acb_vec_is_real(z, g); */
/*     acb_ptr zero; */
/*     arb_ptr w; */
/*     slong a; */

/*     zero = _acb_vec_init(g); */
/*     w = _arb_vec_init(g); */

/*     /\* Set exponentials, etc. *\/ */
/*     acb_theta_ctx_set_z(ctx, zero, 0, prec); */
/*     if (!z_is_zero) */
/*     { */
/*         acb_theta_ctx_set_z(ctx, z, 3, prec); */
/*     } */

/*     /\* Compute distances *\/ */
/*     if (g >= 2) */
/*     { */
/*         for (a = 0; a < n; a++) */
/*         { */
/*             acb_theta_char_get_arb(w, a, g); */
/*             arb_mat_vector_mul_col(w, acb_theta_ctx_cho(ctx), w, lp); */
/*             acb_theta_dist_lat(&acb_theta_ctx_d0(ctx)[a], w, acb_theta_ctx_cho(ctx), lp); */
/*         } */
/*         if (!z_is_real) */
/*         { */
/*             for (a = 0; a < n; a++) */
/*             { */
/*                 acb_theta_char_get_arb(w, a, g); */
/*                 arb_mat_vector_mul_col(w, acb_theta_ctx_cho(ctx), w, lp); */
/*                 _arb_vec_add(w, acb_theta_ctx_vs(ctx) + 3 * g, w, g, lp); */
/*                 acb_theta_dist_lat(&acb_theta_ctx_d(ctx)[a], w, acb_theta_ctx_cho(ctx), lp); */
/*             } */
/*         } */
/*     } */

/*     /\* Set info *\/ */
/*     ctx->t_is_zero = 1; */
/*     ctx->z_is_zero = z_is_zero; */
/*     ctx->z_is_real = z_is_real; */

/*     _acb_vec_clear(zero, g); */
/*     _arb_vec_clear(w, g); */
/* } */
