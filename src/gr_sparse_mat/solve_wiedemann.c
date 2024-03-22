/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "gr_sparse_mat.h"

/* Berlekamp - Massey algorithm */
static int find_min_poly(slong *L, gr_ptr s, slong N, gr_ctx_t ctx)
{
    slong m, n, i, sz;
    slong deg_C, deg_B, deg_T;
    gr_ptr B, C, T;
    gr_ptr c, d_C, d_B;
    int status = GR_SUCCESS;

    sz = ctx->sizeof_elem;
    GR_TMP_INIT3(c, d_C, d_B, ctx);
    GR_TMP_INIT_VEC(B, N, ctx);
    GR_TMP_INIT_VEC(C, N, ctx);
    GR_TMP_INIT_VEC(T, N, ctx);
	
    deg_C = 0, deg_B = 0, deg_T = -1;
    status |= gr_one(d_B, ctx);
    status |= gr_one(B, ctx);
    status |= gr_one(C, ctx);

    *L = 0;
	for (n = 0, m = 1; n < N; n++, m++)
	{
		/* d_C = sum_{i = 0}^L C_i * s_{n-i} */
        status |= gr_set(d_C, GR_ENTRY(s, n, sz), ctx);
		for (i = 1; i <= *L; i++)
			status |= gr_addmul(d_C, GR_ENTRY(C, i, sz), GR_ENTRY(s, n-i, sz), ctx);
        if (gr_is_zero(d_C, ctx) == T_TRUE) continue; /* C and L currently valid */

        /* C(x) = C(x) - (d_C/d_B) x^m B(x); */
        if (*L <= 2*n)
        {
            deg_T = deg_C;
            status |= _gr_vec_set(T, C, deg_C+1, ctx);
        }
        status |= gr_div(c, d_C, d_B, ctx);
        status |= gr_neg(c, c, ctx);
        for (i = 0; i <= deg_B; ++i)
            status |= gr_addmul(GR_ENTRY(C, m+i, sz), GR_ENTRY(B, i, sz), c, ctx);
        
        deg_C = FLINT_MAX(deg_C, deg_B + m);
        while (gr_is_zero(GR_ENTRY(C, deg_C, sz), ctx) == T_TRUE)
            --deg_C;  /* Probably unnecessary */

        if (2 * *L <= n) /* Increase number of errors */
        {
            *L = n + 1 - *L, m = 0;
            status |= gr_set(d_B, d_C, ctx);
            deg_B = deg_T;
            status |= _gr_vec_set(B, T, deg_T+1, ctx);
        }
	}
    /* Reverse C into s */
    for (i = 0; i <= *L; ++i)
        status |= gr_set(GR_ENTRY(s, i, sz),  GR_ENTRY(C, *L-i, sz), ctx);

    GR_TMP_CLEAR_VEC(B, N, ctx);
    GR_TMP_CLEAR_VEC(C, N, ctx);
    GR_TMP_CLEAR_VEC(T, N, ctx);
    GR_TMP_CLEAR3(c, d_C, d_B, ctx);
	return status;
}

/* Compute s_ij=(M^j y)_i for i = 0,...,ns-1, j = 0,...,num-1*/
static int make_sequences(gr_ptr *s, slong ns, slong len, const gr_lil_mat_t M, gr_srcptr b, gr_ctx_t ctx) 
{
    slong i, j, r, sz;
    gr_ptr y, My;
    int status = GR_SUCCESS;

    sz = ctx->sizeof_elem;
    r = gr_sparse_mat_nrows(M, ctx);

    GR_TMP_INIT_VEC(y, r, ctx);
    GR_TMP_INIT_VEC(My, r, ctx);
    status |= _gr_vec_set(y, b, r, ctx);

    for (j = 0; j < len; ++j) 
    {
        if (j > 0)
        {
            status |= gr_lil_mat_mul_vec(My, M, y, ctx);
            status |= _gr_vec_set(y, My, r, ctx);
        }
        for (i = 0; i < ns; ++i)
            status |= gr_set(GR_ENTRY(s[i], j, sz), GR_ENTRY(y, i, sz), ctx);
    }
    GR_TMP_CLEAR_VEC(y, r, ctx);
    GR_TMP_CLEAR_VEC(My, r, ctx);
    return status;
}

/* Compute x = \Sigma_{i = 0}^{L-1} s_i * M^i * b = 0 */
static int make_sum(gr_ptr x, gr_ptr s, slong L, const gr_lil_mat_t M, gr_srcptr b, gr_ctx_t ctx)
{
    slong i, r, sz;
    gr_ptr y, My;
    int status;

    sz = ctx->sizeof_elem;
    r = gr_sparse_mat_nrows(M, ctx);

    GR_TMP_INIT_VEC(y, r, ctx);
    GR_TMP_INIT_VEC(My, r, ctx);
    status = _gr_vec_set(y, b, r, ctx);

    //flint_printf("\t\tScaling\n");
    status = _gr_vec_mul_scalar(x, b, r, s, ctx);
    for (i = 1; i < L; ++i) 
    {
        //flint_printf("\t\tIterating %d\n", i);
        status |= gr_lil_mat_mul_vec(My, M, y, ctx);
        status |= _gr_vec_set(y, My, r, ctx);
        status |= _gr_vec_addmul_scalar(x, y, r, GR_ENTRY(s, i, sz), ctx);
    }
    GR_TMP_CLEAR_VEC(y, r, ctx);
    GR_TMP_CLEAR_VEC(My, r, ctx);
    return status;
}

int gr_lil_mat_solve_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, gr_ctx_t ctx)
{
    slong i, r, c, L, ns, len, sz;
    gr_ptr *s; 
    gr_ptr Mx, coeff;
    int status = GR_SUCCESS;

    // TODO: handle this
    if (x == b)
        return GR_DOMAIN;
    /* TBD: reduce to square */
    if (M->r != M->c)
        return GR_DOMAIN;

    sz = ctx->sizeof_elem;
    r = gr_sparse_mat_nrows(M, ctx);
    c = gr_sparse_mat_ncols(M, ctx);

    if (_gr_vec_is_zero(b, c, ctx) == T_TRUE)
    {
        return _gr_vec_zero(x, c, ctx);
    }

    GR_TMP_INIT_VEC(Mx, r, ctx);
    GR_TMP_INIT(coeff, ctx);

    // Get dimension of sequence to solve
    ns = FLINT_MIN(r, 2);
    len = 2 * r + 1;
    s = flint_malloc(ns * sizeof(gr_ptr));
    for (i = 0; i < ns; ++i)
        GR_TMP_INIT_VEC(s[i], len, ctx);

    //flint_printf("Make sequences\n");
    status |= make_sequences(s, ns, len, M, b, ctx);

    /* Don't have block Berlekamp yet, just try each one */
    for (i = 0; i < ns; ++i)
    {
        /* Get minimal polynomial */
        //flint_printf("Finding minimal polynomial for index %d\n", i);
        status |= find_min_poly(&L, s[i], len, ctx);
        if (gr_is_zero(s[i], ctx) == T_TRUE) continue;

        /* \sum_{j = 0}^L s_ijM^jb = 0 */
        /* => x = -1/s[0]\sum_{j = 0}^{L-1} s_i(j-1) M^jb solves Mx = b */
        //flint_printf("\tMaking sum\n");
        status |= make_sum(x, GR_ENTRY(s[i], 1, sz), L, M, b, ctx);
        //flint_printf("\tInverting\n");
        status |= gr_inv(coeff, s[i], ctx);
        status |= gr_neg(coeff, coeff, ctx);
        //flint_printf("\tScaling\n");
        status |= _gr_vec_mul_scalar(x, x, r, coeff, ctx);

        /* Check if successful */
        //flint_printf("\tChecking result\n");
        status |= gr_lil_mat_mul_vec(Mx, M, x, ctx);
        if (_gr_vec_equal(Mx, b, r, ctx) != T_FALSE)
            break;
    }
    if (i == ns)
        status |= GR_UNABLE;

    GR_TMP_CLEAR_VEC(Mx, r, ctx);
    GR_TMP_CLEAR(coeff, ctx);

    for (i = 0; i < ns; ++i)
        GR_TMP_CLEAR_VEC(s[i], len, ctx);
    flint_free(s);
    return status;
}

int gr_lil_mat_nullvector_wiedemann(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx) 
{
    slong i, L, r, c, ns, len;
    gr_ptr *s; 
    gr_ptr Mx, b;
    int status = GR_SUCCESS;

    // TODO: handle this
     if (M->r != M->c)
        return GR_DOMAIN;

    r = gr_sparse_mat_nrows(M, ctx);
    c = gr_sparse_mat_ncols(M, ctx);

    GR_TMP_INIT_VEC(Mx, r, ctx);
    GR_TMP_INIT_VEC(b, r, ctx);

    // Get dimension of sequence to solve
    ns = FLINT_MIN(r, 2);
    len = 2 * r + 1;
    s = flint_malloc(ns * sizeof(gr_ptr));
    for (i = 0; i < ns; ++i)
        GR_TMP_INIT_VEC(s[i], len, ctx);

    status |= _gr_vec_randtest(x, state, r, ctx);
    status |= gr_lil_mat_mul_vec(b, M, x, ctx);

    status |= make_sequences(s, ns, len, M, b, ctx);

    for (i = 0; i < ns; ++i)
    {
        /* Get minimal polynomial */
        status |= find_min_poly(&L, s[i], len, ctx);

        /* \sum_{j = 0}^L s_ijM^jb = 0 */
        /* => x = \sum_{j = 0}^L s_ijM^jx solves Mx = 0 */
        status |= make_sum(x, s[i], L+1, M, x, ctx);
        status |= gr_lil_mat_mul_vec(Mx, M, x, ctx);
        if 
        (
            _gr_vec_is_zero(x, c, ctx) != T_TRUE || 
            _gr_vec_is_zero(Mx, r, ctx) != T_FALSE
        )
            break;
    }
    if (i == ns)
        status |= GR_UNABLE;

    GR_TMP_CLEAR_VEC(Mx, r, ctx);
    GR_TMP_CLEAR_VEC(b, r, ctx);
    for (i = 0; i < ns; ++i)
        GR_TMP_CLEAR_VEC(s[i], len, ctx);
    flint_free(s);
    return status;
}
