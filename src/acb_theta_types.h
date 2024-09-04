/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_THETA_TYPES_H
#define ACB_THETA_TYPES_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct acb_theta_eld_struct
{
    slong dim, ambient_dim;
    slong * last_coords;
    slong min, mid, max, nr, nl;
    struct acb_theta_eld_struct * rchildren;
    struct acb_theta_eld_struct * lchildren;
    slong nb_pts, nb_border;
    slong * box;
};

typedef struct acb_theta_eld_struct acb_theta_eld_t[1];

struct acb_theta_ctx_tau_struct
{
    acb_mat_struct tau;
    arb_mat_struct Y;
    arb_mat_struct Yinv;
    acb_mat_t exp_tau_div_4;
    acb_mat_t exp_tau_div_2;
    acb_mat_t exp_tau;

    /* g > 1 only */
    arb_mat_struct C;
    arb_mat_struct Cinv;
    acb_mat_t exp_tau_div_4_inv;
    acb_mat_t exp_tau_div_2_inv;
    acb_mat_t exp_tau_inv;
    acb_ptr exp_tau_a_div_2;
    acb_ptr exp_tau_a;
    acb_ptr exp_tau_a_div_2_inv;
    acb_ptr exp_tau_a_inv;
    acb_ptr exp_a_tau_a_div_4;
};

typedef struct acb_theta_ctx_tau_struct acb_theta_ctx_tau_t[1];

typedef struct
{
    slong g;
    acb_ptr z;
    acb_ptr exp_z;
    acb_struct c;
    arb_ptr r;
    arb_struct u;
    arb_struct uinv;
    int is_real;

    /* g > 1 only */
    acb_ptr exp_2z;
    acb_ptr exp_z_inv;
    acb_ptr exp_2z_inv;
    arb_ptr v;
}
acb_theta_ctx_z_struct;

typedef acb_theta_ctx_z_struct acb_theta_ctx_z_t[1];

typedef void (*acb_theta_sum_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, const slong *,
    slong, const acb_t, const slong *, slong, slong, slong, slong);

#ifdef __cplusplus
}
#endif

#endif
