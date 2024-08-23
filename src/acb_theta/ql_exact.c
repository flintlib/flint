/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, */
/*     int sqr, slong prec) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     slong n = 1 << g; */
/*     acb_theta_ctx_t ctx; */
/*     acb_ptr t, rts, zero; */
/*     slong split, nb_steps, nb_ts, nb_rts; */
/*     slong guard = 16; */
/*     int found_t; */

/*     if (nb <= 0) */
/*     { */
/*         return; */
/*     } */

/*     /\* In the worst case, we could use as much as 1 + 5 * nb values of z, */
/*        nb distinct values of t, 4 * n * nb * nb_steps entries of rts, */
/*        and ... values of theta. *\/ */
/*     acb_theta_ctx_init(ctx, 1 + 5 * nb, g); */
/*     t = _acb_vec_init(nb * g); */
/*     zero = _acb_vec_init(g); */

/*     /\* Find out whether to split now. *\/ */
/*     acb_theta_ctx_set_tau(ctx, tau, guard); */
/*     nb_steps = acb_theta_nb_steps(&split, ctx, prec); */

/*     if (split > 0) */
/*     { */
/*         ??? */
/*     } */
/*     else if (nb_steps == 0) */
/*     { */
/*         ??? */
/*     } */

/*     /\* Assume nb_steps is at least 1 and split is 0. *\/ */

/*     nb_rts = 4 * n * nb_steps; /\* not right, we forget the very first step *\/ */
/*     rts = _acb_vec_init(nb * nb_rts); */

/*     /\* Choose vectors t, a common guard, and compute roots. */
/*        In the worst case, there could be as much as nb distinct values of t. *\/ */
/*     t = _acb_vec_init(nb * g); */

/*     found_t = 0; */
/*     for (guard = 16; (guard < prec) && !found_t; guard *= 2) */
/*     { */
/*         /\* Set context at precision guard *\/ */
/*         acb_theta_ctx_set_tau(ctx, tau, guard); */
/*         acb_theta_ctx_set_z(ctx, zero, 0, guard); */
/*         /\* For each z, we will possibly store t, 2t, z, z+t, z+2t *\/ */
/*         for (j = 0; j < nb; j++) */
/*         { */
/*             acb_theta_ctx_set_z(ctx, z, 3 + 5 * j, guard); */
/*         } */

/*         _acb_vec_zero(t, nb * g); */
/*         nb_ts = 1; /\* Always have the zero vector. *\/ */
/*         found_t = 1; */
/*         for (j = 0; (j < nb) && found_t; j++) */
/*         { */
/*             found_t = 0; */
/*             for (k = 0; (k < nb_ts) && !found_t; k++) */
/*             { */
/*                 /\* Attempt to compute roots for zj, tk at current guard. *\/ */
/*                 /\* Set z + t, z + 2t; we already know 0, t, 2t. *\/ */
/*                 /\* Nothing to be done if k = 0, ie t corresponds to zero vector. *\/ */
/*                 acb_theta_ctx_set_ql_translates(ctx, j, k, guard); */
/*                 found_t = acb_theta_ql_roots_z(rts + j * nb_rts, ctx, j, k, guard); */
/*             } */
/*             if (found_t) /\* We are happy. k was increased before exiting loop. *\/ */
/*             { */
/*                 ctx->t_indices[j] = k - 1; */
/*             } */
/*             else */
/*             { */
/*                 /\* We did not find any t that works. Pick a few new ones *\/ */
/*                 /\* k is nb_ts *\/ */
/*                 for (l = 0; (l < 10) && !found_t; l++) */
/*                 { */
/*                     for (m = 0; m < g; m++) */
/*                     { */
/*                         arb_urandom(acb_realref(&t[g * k + k])); */
/*                     } */
/*                     /\* Use slots 1 + 5 * k, 2 + 5 * k for t and 2t *\/ */
/*                     acb_theta_ctx_set_t(ctx, t, 1 + 5 * k, guard); */
/*                     found_t = acb_theta_ql_roots_0(ctx, k, guard); */
/*                     if (!found_t) */
/*                     { */
/*                         continue; */
/*                     } */
/*                     /\* This choice of t just worked for roots at t, 2t. *\/ */
/*                     acb_theta_ctx_set_ql_translates(ctx, j, k, guard); */
/*                     found_t = acb_theta_ql_roots_z(rts + j * nb_rts, ctx, j, k, guard); */
/*                 } */
/*                 if (found_t) */
/*                 { */
/*                     ctx->t_indices[j] = k; */
/*                 } */
/*                 nb_ts += 1; */
/*             } */
/*             /\* Now either we found t, and we happily continue the loop on j; */
/*                or we didn't, exit the loop on j, and increase guard. *\/ */
/*         } */
/*     } */

/*     if (found_t) */
/*     { */
/*         /\* We have our choice of ts, guards, and all the relevant roots have been computed. *\/ */

/*         /\* Strip z, tau of error bounds *\/ */

/*         /\* Initialize theta high up. *\/ */
/*         acb_theta_ql_initialize(th, ...); */

/*         /\* Duplication steps *\/ */
/*         /\* Precision may decrease slightly at each step, but doesn't matter too much. *\/ */

/*         /\* Final duplication step is special. *\/ */

/*         /\* Finally, add error bounds coming from radii of tau, z, nothing fancy. *\/ */
/*     } */
/* } */
