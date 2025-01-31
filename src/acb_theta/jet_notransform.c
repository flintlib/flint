/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

/* We use the formula:
   theta_00(z, tau) = sum_a theta_{a,0}(2z, 4tau) */

/* static void */
/* acb_theta_jet_00_notransform(acb_ptr th, acb_srcptr zs, slong nb, */
/*     const acb_mat_t tau, slong ord, slong prec) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     slong n = 1 << g; */
/*     slong nbder = acb_theta_jet_nb(ord, g); */
/*     acb_ptr new_zs, aux; */
/*     acb_mat_t new_tau; */

/*     new_zs = _acb_vec_init(nb * g); */
/*     aux = _acb_vec_init(nb * n * nbder); */
/*     acb_mat_init(new_tau, g, g); */

/*     _acb_vec_scalar_mul_2exp_si(new_zs, zs, nb * g, 1); */
/*     acb_mat_scalar_mul_2exp_si(new_tau, tau, 2); */

/*     acb_theta_jet_a0_notransform(aux, new_zs, nb, new_tau, ord, prec); */
/*     _acb_vec_zero(th, nb * nbder); */
/*     for (j = 0; j < nb; j++) */
/*     { */
/*         for (k = 0; k < nbder; k++) */
/*         { */
/*             for (a = 0; a < n; a++) */
/*             { */
/*                 acb_add(&th[j * nbder + k], &th[j * nbder + k], */
/*                     &aux[j * n * nbder + a * nbder + k], prec); */
/*                 /\* Need to multiply by power of 2 *\/ */
/*             } */
/*         } */
/*     } */

/*     _acb_vec_clear(new_zs, nb * g); */
/*     _acb_vec_clear(aux, nb * n * nbder); */
/*     acb_mat_clear(new_tau); */
/* } */

/* /\* We use the formula: */
/*    theta_ab(z, tau) = exp(pi i a^T tau a/4) exp(2 pi i a^T (z + b/2)) */
/*            theta_00(z + tau a/2 + b/2, tau) *\/ */

/* static void */
/* acb_theta_jet_one_notransform(acb_ptr th, acb_srcptr zs, slong nb, */
/*     const acb_mat_t tau, slong ord, ulong ab, slong prec) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     slong nbth = acb_theta_jet_nb(ord, g); */
/*     ulong b = ab % (1 << g); */
/*     ulong a = ab >> g; */
/*     acb_ptr new_zs, v, w, aux; */
/*     arb_ptr u; */
/*     acb_t c, x; */
/*     slong j; */

/*     new_zs = _acb_vec_init(nb * g); */
/*     v = _acb_vec_init(g); */
/*     w = _acb_vec_init(g); */
/*     aux = _acb_vec_init(nbth); */
/*     u = _arb_vec_init(g); */
/*     acb_init(c); */
/*     acb_init(x); */

/*     acb_theta_char_get_acb(v, a, g); */
/*     acb_mat_vector_mul_col(v, tau, v, prec); /\* tau.a/2 *\/ */
/*     acb_theta_char_get_acb(w, b, g); */
/*     _acb_vec_add(w, v, w, g, prec); */
/*     for (j = 0; j < nb; j++) */
/*     { */
/*         _acb_vec_add(new_zs + j * g, zs + j * g, w, g, prec); */
/*     } */

/*     acb_theta_jet_00_notransform(th, new_zs, nb, tau, ord, prec); */

/*     acb_theta_char_dot_acb(c, a, v, g, prec); */
/*     for (j = 0; j < nb; j++) */
/*     { */
/*         acb_theta_char_get_acb(w, b, g); */
/*         _acb_vec_add(w, w, zs + j * g, g, prec); */
/*         acb_theta_char_dot_acb(x, a, w, g, prec); */
/*         acb_mul_2exp_si(x, x, 1); */
/*         acb_add(x, x, c, prec); */
/*         acb_exp_pi_i(x, x, prec); */
/*         _acb_vec_scalar_mul(th + j * nbth, th + j * nbth, nbth, x, prec); */
/*     } */

/*     if (ord > 0) */
/*     { */
/*         acb_theta_char_get_arb(u, a, g); */
/*         _arb_vec_scalar_mul_2exp_si(u, u, g, 1); */
/*         acb_theta_jet_exp_pi_i(aux, u, ord, g, prec); */
/*         for (j = 0; j < nb; j++) */
/*         { */
/*             acb_theta_jet_mul(th + j * nbth, th + j * nbth, aux, ord, g, prec); */
/*         } */
/*     } */

/*     _acb_vec_clear(new_zs, nb * g); */
/*     _acb_vec_clear(v, g); */
/*     _acb_vec_clear(w, g); */
/*     _acb_vec_clear(aux, nbth); */
/*     _arb_vec_clear(u, g); */
/*     acb_clear(c); */
/*     acb_clear(x); */
/* } */

/* void */
/* acb_theta_jet_notransform(acb_ptr th, acb_srcptr zs, slong nb, */
/*     const acb_mat_t tau, slong ord, ulong ab, int all, slong prec) */
/* { */
    
/* } */
