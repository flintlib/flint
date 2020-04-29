/*
    Copyright (C) 2010 Fredrik Johansson

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
#include "nmod_sparse_mat.h"

/* Berlekamp - Massey algorithm */
static slong find_min_poly(mp_limb_t *s, slong N, nmod_t mod)
{
    slong L = 0, m, n, i;
    mp_limb_t c, d_C, d_B = 1;
	
    slong deg_C = 0, deg_B = 0, deg_T = -1;
    mp_limb_t *B, *C, *T;
    B = flint_calloc(N, sizeof(*B));
    C = flint_calloc(N, sizeof(*C));
    T = flint_calloc(N, sizeof(*T));
    B[0] = C[0] = UWORD(1);

	for (n = 0, m = 1; n < N; n++, m++)
	{
		/* d_C = sum_{i = 0}^L C_i * s_{n-i} */
		d_C = s[n];
		for (i = 1; i <= L; i++)
			d_C = nmod_addmul(d_C, C[i], s[n-i], mod);
        if (d_C == 0) continue; /* C and L currently valid */

        /* C(x) = C(x) - (d_C/d_B) x^m B(x); */
        if (L <= 2*n) deg_T = deg_C, memcpy(T, C, (deg_C+1)*sizeof(*T)); /* T(X) = C(X) */
        c = nmod_neg(nmod_div(d_C, d_B, mod), mod);
        for (i = 0; i <= deg_B; ++i)
            C[m+i] = nmod_addmul(C[m+i], B[i], c, mod);
        deg_C = FLINT_MAX(deg_C, deg_B + m);
        while (C[deg_C] == UWORD(0)) --deg_C;  /* Probably unnecessary */

        if (2*L <= n) /* Increase number of errors */
        {
            L = n + 1 - L, m = 0;
            d_B = d_C, deg_B = deg_T;
            memcpy(B, T, (deg_T+1)*sizeof(*B)); /* B(x) = C(x) */
        }
	}
    /* Reverse C into s */
    for (i = 0; i <= L; ++i) s[i] = C[L-i];
    flint_free(B);
    flint_free(C);
    flint_free(T);
	return L;
}

/* Compute s_ij=(M^j y)_i for i = 0,...,ns-1, j = 0,...,num-1*/
static void make_sequences(mp_limb_t **s, slong ns, slong len, const nmod_sparse_mat_t M, mp_srcptr b) 
{
    slong i, j;
    mp_ptr y, My;
    y = _nmod_vec_init(M->r);
    My = _nmod_vec_init(M->r);
    memcpy(y, b, M->r*sizeof(*y));
    for (j = 0; j < len; ++j) 
    {
        if (j > 0) nmod_sparse_mat_mul_vec(My, M, y), memcpy(y, My, M->r*sizeof(*y));
        for (i = 0; i < ns; ++i) s[i][j] = y[i];
    }
    _nmod_vec_clear(y);
    _nmod_vec_clear(My);
}

/* Compute x = \Sigma_{i = 0}^{L-1} s_i * M^i * b = 0 */
static void make_sum(mp_ptr x, mp_limb_t *s, slong L, const nmod_sparse_mat_t M, mp_srcptr b)
{
    slong i;
    mp_ptr y, My;
    y = _nmod_vec_init(M->r);
    My = _nmod_vec_init(M->r);
    memcpy(y, b, M->r*sizeof(*y));
    _nmod_vec_scalar_mul_nmod(x, b, M->r, s[0], M->mod);
    for (i = 1; i < L; ++i) 
    {
        nmod_sparse_mat_mul_vec(My, M, y), memcpy(y, My, M->r*sizeof(*y));
        _nmod_vec_scalar_addmul_nmod(x, y, M->r, s[i], M->mod);
    }
    _nmod_vec_clear(y);
    _nmod_vec_clear(My);
}

int nmod_sparse_mat_solve_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b)
{
    slong i, L, ret = 0, ns = FLINT_MIN(M->r, 2), len = 2*M->r + 1;
    mp_limb_t **s; 
    mp_ptr Mx;
    if (M->r != M->c) return 0; /* TBD: reduce to square */
    if (_nmod_vec_is_zero(b, M->c))
    {
        _nmod_vec_zero(x, M->c);
        return 1;
    }

    Mx = _nmod_vec_init(M->r);
    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = flint_malloc(len*sizeof(*s[i]));
    
    make_sequences(s, ns, len, M, b);

    /* Don't have block Berlekamp yet, just try each one */
    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, M->mod);
        if (s[i][0]==0) continue;

        /* If \sum_{j = 0}^L s_ijM^jb = 0 => x = -1/s[0]\sum_{j = 0}^{L-1} s_i(j-1) M^jb solves Mx = b */
        make_sum(x, s[i]+1, L, M, b);
        _nmod_vec_scalar_mul_nmod(x, x, M->r, nmod_neg(nmod_inv(s[i][0], M->mod), M->mod), M->mod);

        /* Check if successful */
        nmod_sparse_mat_mul_vec(Mx, M, x);
        ret = _nmod_vec_equal(Mx, b, M->r);
    }

    _nmod_vec_clear(Mx);
    for (i = 0; i < ns; ++i) flint_free(s[i]);
    flint_free(s);
    return ret;
}

int nmod_sparse_mat_nullvector_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state) 
{
    slong i, L, ret = 0, ns = FLINT_MIN(M->r, 2), len = 2*M->r + 1;
    mp_limb_t **s; 
    mp_ptr Mx, b;
    Mx = _nmod_vec_init(M->r);
    b = _nmod_vec_init(M->r);

    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = flint_malloc(len*sizeof(*s[i]));

    _nmod_vec_randtest(x, state, M->r, M->mod);
    nmod_sparse_mat_mul_vec(b, M, x);

    if (M->r != M->c) return 0; /* TBD: reduce to square */
    make_sequences(s, ns, len, M, b);

    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, M->mod);

        /* \sum_{j = 0}^L s_ijM^jb = 0 => x = \sum_{j = 0}^L s_ijM^jx solves Mx = 0 */
        make_sum(x, s[i], L+1, M, x);
        nmod_sparse_mat_mul_vec(Mx, M, x);
        ret = !_nmod_vec_is_zero(x, M->c) && _nmod_vec_is_zero(Mx, M->r);
    }

    _nmod_vec_clear(Mx);
    _nmod_vec_clear(b);
    for (i = 0; i < ns; ++i) flint_free(s[i]);
    flint_free(s);
    return ret;
}