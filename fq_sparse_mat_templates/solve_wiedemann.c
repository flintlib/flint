/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

/* Berlekamp - Massey algorithm */
static slong find_min_poly(TEMPLATE(T, struct) *s, slong N, const TEMPLATE(T, ctx_t) ctx)
{
    slong L = 0, m, n, i;
    slong deg_C = 0, deg_B = 0, deg_T = -1;
    TEMPLATE(T, t) cinv, cc, d_C, d_B;
    TEMPLATE(T, struct) *B, *C, *T;
	
    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, init) (d_C, ctx);
    TEMPLATE(T, init) (d_B, ctx);
    B = _TEMPLATE(T, vec_init) (N, ctx);
    C = _TEMPLATE(T, vec_init) (N, ctx);
    T = _TEMPLATE(T, vec_init) (N, ctx);
    TEMPLATE(T, one) (&B[0], ctx);
    TEMPLATE(T, one) (&C[0], ctx);
    TEMPLATE(T, one) (d_B, ctx);

	for (n = 0, m = 1; n < N; n++, m++)
	{
		/* d_C = sum_{i = 0}^L C_i * s_{n-i} */
		TEMPLATE(T, set) (d_C, &s[n], ctx);
		for (i = 1; i <= L; i++)
        {
            TEMPLATE(T, mul) (cc, &C[i], &s[n-i], ctx);
			TEMPLATE(T, add) (d_C, d_C, cc, ctx);
        }
        if (TEMPLATE(T, is_zero)(d_C, ctx)) continue; /* C and L currently valid */

        /* C(x) = C(x) - (d_C/d_B) x^m B(x); */
        if (L <= 2*n) deg_T = deg_C, _TEMPLATE(T, vec_set) (T, C, deg_C+1, ctx); /* T(X) = C(X) */
        TEMPLATE(T, div) (cinv, d_C, d_B, ctx);
        TEMPLATE(T, neg) (cinv, cinv, ctx);
        for (i = 0; i <= deg_B; ++i)
        {
            TEMPLATE(T, mul) (cc, &B[i], cinv, ctx);
            TEMPLATE(T, add) (&C[m+i], &C[m+i], cc, ctx);
        }
        deg_C = FLINT_MAX(deg_C, deg_B + m);
        while (TEMPLATE(T, is_zero) (&C[deg_C], ctx)) --deg_C;  /* Probably unnecessary */

        if (2*L <= n) /* Increase number of errors */
        {
            L = n + 1 - L, m = 0;
            TEMPLATE(T, set) (d_B, d_C, ctx), deg_B = deg_T;
            _TEMPLATE(T, vec_set) (B, T, deg_T+1, ctx);
        }
	}
    /* Reverse C into s */
    for (i = 0; i <= L; ++i) TEMPLATE(T, set) (&s[i], &C[L-i], ctx);

    TEMPLATE(T, clear) (cinv, ctx);
    TEMPLATE(T, clear) (cc, ctx);
    TEMPLATE(T, clear) (d_C, ctx);
    TEMPLATE(T, clear) (d_B, ctx);
    _TEMPLATE(T, vec_clear) (B, N, ctx);
    _TEMPLATE(T, vec_clear) (C, N, ctx);
    _TEMPLATE(T, vec_clear) (T, N, ctx);
	return L;
}

/* Compute s_ij=(A^j y)_i for i = 0,...,ns-1, j = 0,...,num-1*/
static void make_sequences(TEMPLATE(T, struct) **s, slong ns, slong len, const TEMPLATE(T, sparse_mat_t) A, TEMPLATE(T, struct) **y, const TEMPLATE(T, ctx_t) ctx) 
{
    slong iter, i, j;
    for (i = iter = 0; iter < len; ++iter, i = 1-i) 
    {
        if (iter > 0) TEMPLATE(T, sparse_mat_mul_vec) (y[i], A, y[1-i], ctx);
        for (j = 0; j < ns; ++j) TEMPLATE(T, set) (&s[j][iter], &y[i][j], ctx);
    }
}

/* Compute x = \Sigma_{i = 0}^{L-1} s_i * A^i * b = 0 */
static void make_sum(TEMPLATE(T, struct) *x, TEMPLATE(T, struct) *s, slong L, const TEMPLATE(T, sparse_mat_t) A, TEMPLATE(T, struct) **y, const TEMPLATE(T, ctx_t) ctx)
{
    slong iter, i;
    _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (x, y[0], A->r, &s[0], ctx);
    for (i = iter = 1; iter < L; ++iter, i = 1-i) 
    {
        TEMPLATE(T, sparse_mat_mul_vec) (y[i], A, y[1-i], ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (x, y[i], A->r, &s[iter], ctx);
    }
}

int TEMPLATE(T, sparse_mat_solve_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) A, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, L, ret = 0, ns = FLINT_MIN(A->r, 2), len = 2*A->r + 1;
    TEMPLATE(T, t) cc;  
    TEMPLATE(T, struct) **s, *y[2];
    if (A->r != A->c) return 0; /* TBD: reduce to square */
    if (_TEMPLATE(T, vec_is_zero) (b, A->r, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, A->c, ctx);
        return 1;
    }
    TEMPLATE(T, init) (cc, ctx);
    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = _TEMPLATE(T, vec_init) (len, ctx);
    for (i = 0; i < 2; ++i) y[i] = _TEMPLATE(T, vec_init) (A->r, ctx);
    
    _TEMPLATE(T, vec_set) (y[0], b, A->r, ctx);

    /* Make some number of sequences to be tried */
    make_sequences(s, ns, len, A, y, ctx);

    /* Try both sequences */
    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, ctx);
        if (TEMPLATE(T, is_zero) (&s[i][0], ctx)) continue;

        /* If \sum_{j = 0}^L s_ijA^jb = 0 => x = -1/s[0]\sum_{j = 0}^{L-1} s_i(j-1) A^jb solves Ax = b */
        _TEMPLATE(T, vec_set) (y[0], b, A->r, ctx);
        make_sum(x, s[i]+1, L, A, y, ctx);
        TEMPLATE(T, inv) (cc, &s[i][0], ctx);
        TEMPLATE(T, neg) (cc, cc, ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (x, x, A->r, cc, ctx);

        /* Check if successful */
        TEMPLATE(T, sparse_mat_mul_vec) (y[0], A, x, ctx);
        ret = _TEMPLATE(T, vec_equal) (y[0], b, A->r, ctx);
    }

    TEMPLATE(T, clear) (cc, ctx);
    for (i = 0; i < ns; ++i) _TEMPLATE(T, vec_clear) (s[i], len, ctx);
    for (i = 0; i < 2; ++i) _TEMPLATE(T, vec_clear) (y[i], A->r, ctx);
    flint_free(s);
    return ret;
}

int TEMPLATE(T, sparse_mat_nullvector_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) A, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i, L, ret = 0, ns = FLINT_MIN(A->r, 2), len = 2*A->r + 1;
    TEMPLATE(T, struct) **s, *y[3]; 

    if (A->r != A->c) return 0; /* TBD: reduce to square */

    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = _TEMPLATE(T, vec_init) (len, ctx);
    for (i = 0; i < 3; ++i) y[i] = _TEMPLATE(T, vec_init) (A->r, ctx);

    _TEMPLATE(T, vec_randtest) (y[0], state, A->r, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (y[1], A, y[0], ctx);

    make_sequences(s, ns, len, A, &y[1], ctx);

    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, ctx);

        /* \sum_{j = 0}^L s_ijA^jb = 0 => x = \sum_{j = 0}^L s_ijA^jx solves Ax = 0 */
        _TEMPLATE(T, vec_set) (y[1], y[0], A->r, ctx);
        make_sum(x, s[i], L+1, A, &y[1], ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (y[1], A, x, ctx);
        ret = !_TEMPLATE(T, vec_is_zero) (x, A->c, ctx) && _TEMPLATE(T, vec_is_zero) (y[1], A->r, ctx);
    }
    for (i = 0; i < ns; ++i) _TEMPLATE(T, vec_clear) (s[i], len, ctx);
    for (i = 0; i < 3; ++i) _TEMPLATE(T, vec_clear) (y[i], A->r, ctx);
    flint_free(s);
    return ret;
}

#endif
