/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_lil_mat_solve_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, flint_rand_t state, gr_ctx_t ctx)
{
    slong j, r, c, iter;
    gr_lil_mat_t Mt;
    gr_ptr v[2], Mv, Av, Mtb;
    gr_ptr vtAv[2], AvtAv, vMtb, tmp;
    int status = GR_SUCCESS;

    // TODO: handle this
    if (x == b)
        return GR_DOMAIN;

    r = gr_sparse_mat_nrows(M, ctx);
    c = gr_sparse_mat_ncols(M, ctx);

    // flint_printf("M = "); gr_lil_mat_print_nz(M, ctx); flint_printf("\n");
    status |= _gr_vec_zero(x, c, ctx);
    if (_gr_vec_is_zero(b, c, ctx) == T_TRUE)
        return GR_SUCCESS; // Trivial solution works

    /* We assume that M is not symmetric, and work with A = M^t M */
    gr_lil_mat_init(Mt, c, r, ctx);
    status |= gr_lil_mat_transpose(Mt, M, ctx);
    //flint_printf("M^T = "); gr_lil_mat_print_nz(Mt, ctx); flint_printf("\n");

    /* Construct auxiliary values and vectors */
    /* Rather than storing the whole sequence of values v_j, we alternate between two vectors */
    GR_TMP_INIT5(vtAv[0], vtAv[1], AvtAv, vMtb, tmp, ctx);
    GR_TMP_INIT_VEC(v[0], c, ctx);
    GR_TMP_INIT_VEC(v[1], c, ctx);
    GR_TMP_INIT_VEC(Mv, r, ctx);
    GR_TMP_INIT_VEC(Av, c, ctx);
    GR_TMP_INIT_VEC(Mtb, c, ctx);
    
    /* Make 0th vector random (and -1st vector trivial) */
    /*_gr_vec_set(v[0], Mtb, M->c);
    for (j = 0; j < M->c; ++j) v[0][j] = n_randint(state, M->mod.n); */
    status |= gr_lil_mat_mul_vec(Mtb, Mt, b, ctx);
    // flint_printf("\tM^Tb = "); _gr_vec_print(Mtb, c, ctx); flint_printf("\n");
    status |= _gr_vec_randtest(v[0], state, c, ctx);
    status |= _gr_vec_zero(v[1], c, ctx);
    status |= gr_one(vtAv[1], ctx);
    for (iter = j = 0; ; j = 1-j, ++iter)
    {
        // flint_printf("\n\niter = %d\n", iter);
        /* Compute A v_j and check if it is orthogonal to v_j */
        // flint_printf("\n\tv[%d] = ", j); _gr_vec_print(v[j], c, ctx); flint_printf("\n");
        status |= gr_lil_mat_mul_vec(Mv, M, v[j], ctx);
        // flint_printf("\tMv[%d] = ", j); _gr_vec_print(Mv, r, ctx); flint_printf("\n");
        status |= gr_lil_mat_mul_vec(Av, Mt, Mv, ctx);
        // flint_printf("\tAv[%d] = ", j); _gr_vec_print(Av, c, ctx); flint_printf("\n");
        status |= _gr_vec_dot(vtAv[j], NULL, 0, v[j], Av, c, ctx);
        // flint_printf("\tv[%d]^TAv[%d] = ", j, j); gr_println(vtAv[j], ctx);
        if (gr_is_zero(vtAv[j], ctx) == T_TRUE) break; /* Can't make any more progress */

        /* Update putative solution by (<v_j, M^T b>/<v_j, Av_j>) * v_j */
        status |= _gr_vec_dot(vMtb, NULL, 0, v[j], Mtb, c, ctx);
        // flint_printf("\n\tv[%d]M^Tb = ", j); gr_println(vMtb, ctx);
        status |= gr_div(vMtb, vMtb, vtAv[j], ctx);
        // flint_printf("\tv[%d]M^Tb/(v[%d]^TAv[%d]) = ", j, j, j); gr_println(vMtb, ctx);
        status |= _gr_vec_addmul_scalar(x, v[j], c, vMtb, ctx);
        // flint_printf("\tx = x + ([%d]M^Tb/(v[%d]^TAv[%d])) v[%d] = ", j, j, j); _gr_vec_print(x, c, ctx); flint_printf("\n");

        /* v_{j+1} = MtMv - alpha*v_j - beta*v_{j-1}, where */
        /*    alpha = <Mv_j, Mv_j>/delta_j, and */
        /*    beta = delta_j/delta_{j-1} */
        status |= _gr_vec_dot(AvtAv, NULL, 0, Av, Av, c, ctx);
        // flint_printf("\t(Av[%d])^TAv[%d] = ", j, j); gr_println(AvtAv, ctx);
        status |= gr_div(tmp, vtAv[j], vtAv[1-j], ctx);
        // flint_printf("\t((Av[%d])^TAv[%d])/((Av[%d])^TAv[%d]) = ", j, j, 1-j, 1-j); gr_println(tmp, ctx);
        status |= gr_neg(tmp, tmp, ctx);
        // flint_printf("\t-((Av[%d])^TAv[%d])/((Av[%d])^TAv[%d]) = ", j, j, 1-j, 1-j); gr_println(tmp, ctx);
        status |= _gr_vec_mul_scalar(v[1-j], v[1-j], M->c, tmp, ctx);
        // flint_printf("\tv[%d] = -((Av[%d])^TAv[%d])/((Av[%d])^TAv[%d]) v[%d] = ", 1-j, j, j, 1-j, 1-j, 1-j); _gr_vec_print(v[1-j], c, ctx); flint_printf("\n");
        
        status |= gr_div(tmp, AvtAv, vtAv[j], ctx);
        // flint_printf("\n\t((Av[%d])^TAv[%d])/(v[%d]^TAv[%d]) = ", j, j, j, j); gr_println(tmp, ctx);
        status |= gr_neg(tmp, tmp, ctx);
        // flint_printf("\t-((Av[%d])^TAv[%d])/(v[%d]^TAv[%d]) = ", j, j, j, j); gr_println(tmp, ctx);
        status |= _gr_vec_addmul_scalar(v[1-j], v[j], M->c, tmp, ctx);
        // flint_printf("\tv[%d] = v[%d] - ((Av[%d])^TAv[%d])/(v[%d]^TAv[%d]) v[%d] = ", 1-j, 1-j, j, j, j, j, j); _gr_vec_print(v[1-j], c, ctx); flint_printf("\n");
        status |= _gr_vec_add(v[1-j], v[1-j], Av, M->c, ctx);
        // flint_printf("\tv[%d] = v[%d] + Av[%d] = ", 1-j, 1-j, j); _gr_vec_print(v[1-j], c, ctx); flint_printf("\n");
    }
    /* Check result */
    // flint_printf("x = ", j); _gr_vec_print(x, c, ctx); flint_printf("\n");
    status |= gr_lil_mat_mul_vec(Mv, M, x, ctx);
    // flint_printf("Mx = ", j); _gr_vec_print(Mv, r, ctx); flint_printf("\n");
    status |= gr_lil_mat_mul_vec(Av, Mt, Mv, ctx);
    // flint_printf("Ax = ", j); _gr_vec_print(Av, r, ctx); flint_printf("\n");
    if (_gr_vec_equal(Av, Mtb, c, ctx) == T_FALSE)
        status = GR_UNABLE;

    /* Clear auxiliary vectors and transpose */
    gr_lil_mat_clear(Mt, ctx);
    GR_TMP_CLEAR5(vtAv[0], vtAv[1], AvtAv, vMtb, tmp, ctx);
    GR_TMP_CLEAR_VEC(v[0], c, ctx);
    GR_TMP_CLEAR_VEC(v[1], c, ctx);
    GR_TMP_CLEAR_VEC(Mv, r, ctx);
    GR_TMP_CLEAR_VEC(Av, c, ctx);
    GR_TMP_CLEAR_VEC(Mtb, c, ctx);
    return status;
}

int gr_lil_mat_nullvector_lanczos(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx) 
{
    int status = GR_SUCCESS;
    gr_ptr x2, b;
    GR_TMP_INIT_VEC(x2, M->c, ctx);
    GR_TMP_INIT_VEC(b, M->r, ctx);

    status |= _gr_vec_randtest(x, state, M->c, ctx);
    status |= gr_lil_mat_mul_vec(b, M, x, ctx);
    status |= gr_lil_mat_solve_lanczos(x2, M, b, state, ctx);

    if (status == GR_SUCCESS) 
    {
        status |= _gr_vec_sub(x, x, x2, M->c, ctx);
        if (_gr_vec_is_zero(x, M->c, ctx) == T_TRUE)
            status = GR_TEST_FAIL;
        else
        {
            status |= gr_lil_mat_mul_vec(b, M, x, ctx);
            if (_gr_vec_is_zero(b, M->r, ctx) == T_FALSE)
                status = GR_DOMAIN;
        }
    }
    GR_TMP_CLEAR_VEC(x2, M->c, ctx);
    GR_TMP_CLEAR_VEC(b, M->r, ctx);
    return status;
}
