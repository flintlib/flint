/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "nfloat.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
fmpz_lll_is_reduced_mpfr(const fmpz_mat_t B, const fmpz_lll_t fl,
                         flint_bitcnt_t prec)
{
    gr_mat_t A, R, V, Wu, Wd, bound;
    gr_ptr du, dd;
    gr_ptr s, norm, ti, tj, tmp;
    gr_ptr zero, one;
    slong n;
    int status = GR_SUCCESS;

    gr_ctx_t ctx, ctx_rnd_floor, ctx_rnd_ceil;

#define ENTRY(vv,ii) GR_ENTRY(vv,ii,ctx->sizeof_elem)
#define MAT_ENTRY(mm,ii,jj) GR_MAT_ENTRY(mm,ii,jj,ctx->sizeof_elem)

    /* Technical detail: we currently require prec >= 53 so that
       conversion from double is exact, because there is no atomic
       gr_mul_d or gr_cmp_d respecting directed rounding. */
    prec = FLINT_MAX(prec, 64);

    status |= nfloat_ctx_init(ctx, prec, 0);
    status |= nfloat_ctx_init(ctx_rnd_floor, prec, NFLOAT_RND_FLOOR);
    status |= nfloat_ctx_init(ctx_rnd_ceil, prec, NFLOAT_RND_CEIL);

    if (status != GR_SUCCESS)
        return 0;

    GR_TMP_INIT2(zero, one, ctx);
    status |= gr_one(one, ctx);

    if (fl->rt == Z_BASIS)
    {
        /* NOTE: this algorithm should *not* be changed */
        slong i, j, k, m;
        gr_mat_t Q, bound2, bound3, boundt, mm, rm, mn, rn, absR;

        if (B->r == 0 || B->r == 1)
        {
            GR_TMP_CLEAR2(zero, one, ctx);
            return 1;
        }

        m = B->c;
        n = B->r;

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(Q, n, m, ctx);
        gr_mat_init(R, n, n, ctx);
        gr_mat_init(V, n, n, ctx);

        GR_TMP_INIT5(s, norm, ti, tj, tmp, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                status |= gr_set_fmpz(MAT_ENTRY(A, j, i), fmpz_mat_entry(B, i, j), ctx);
            }
        }

        /* Here Q is transposed wrt the original double code to allow using vector
           operations. To do: use the same optimization in the double code. */
        for (k = 0; k < n; k++)
        {
            for (j = 0; j < m; j++)
            {
                status |= gr_set(MAT_ENTRY(Q, k, j), MAT_ENTRY(A, j, k), ctx);
            }
            for (i = 0; i < k; i++)
            {
                status |= _gr_vec_dot(s, NULL, 0, MAT_ENTRY(Q, i, 0), MAT_ENTRY(Q, k, 0), m, ctx);
                status |= gr_set(MAT_ENTRY(R, i, k), s, ctx);
                status |= _gr_vec_submul_scalar(MAT_ENTRY(Q, k, 0), MAT_ENTRY(Q, i, 0), m, s, ctx);
            }

            /* todo: optimized vec_norm2 (or optimize self-dot for nfloat, arf) */
            status |= _gr_vec_dot(s, NULL, 0, MAT_ENTRY(Q, k, 0), MAT_ENTRY(Q, k, 0), m, ctx);
            status |= gr_sqrt(s, s, ctx);
            status |= gr_set(MAT_ENTRY(R, k, k), s, ctx);

            if (gr_is_zero(s, ctx) == T_FALSE)
            {
                status |= gr_inv(s, s, ctx);
                status |= _gr_vec_mul_scalar(MAT_ENTRY(Q, k, 0), MAT_ENTRY(Q, k, 0), m, s, ctx);
            }
        }

        gr_mat_clear(Q, ctx);

        for (j = n - 1; j >= 0; j--)
        {
            status |= gr_inv(MAT_ENTRY(V, j, j), MAT_ENTRY(R, j, j), ctx);
            for (i = j + 1; i < n; i++)
            {
                status |= _gr_vec_dot(MAT_ENTRY(V, i, j), MAT_ENTRY(V, i, j), 0,
                    MAT_ENTRY(V, i, j + 1), MAT_ENTRY(R, j, j + 1), n - (j + 1), ctx);
                status |= gr_neg(MAT_ENTRY(V, i, j), MAT_ENTRY(V, i, j), ctx);
                status |= gr_mul(MAT_ENTRY(V, i, j), MAT_ENTRY(V, i, j), MAT_ENTRY(V, j, j), ctx);
            }
        }
        status |= gr_mat_transpose(V, V, ctx);

        gr_mat_init(Wu, n, n, ctx);
        gr_mat_init(Wd, n, n, ctx);
        du = gr_heap_init_vec(n, ctx);
        dd = gr_heap_init_vec(n, ctx);

        status |= gr_mat_mul(Wd, R, V, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(dd, i), MAT_ENTRY(Wd, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(Wu, R, V, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(du, i), MAT_ENTRY(Wu, i, i), 1, ctx_rnd_ceil);
        }
        status |= gr_zero(norm, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_zero(s, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(s, s, tmp, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_abs(ti, ENTRY(dd, i), ctx_rnd_ceil);
                    status |= gr_abs(tj, ENTRY(du, i), ctx_rnd_ceil);
                    status |= gr_max(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(s, s, tmp, ctx_rnd_ceil);
                }
            }
            status |= gr_max(norm, norm, s, ctx_rnd_ceil);
        }

        if (status != GR_SUCCESS || gr_lt(norm, one, ctx) != T_TRUE)
            goto fail_clear_all;

        gr_mat_init(bound, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(dd, i), MAT_ENTRY(Wd, i, i), 2, ctx_rnd_floor);
        }
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(du, i), MAT_ENTRY(Wu, i, i), 2, ctx_rnd_ceil);
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j > i)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
                else if (j < i)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_abs(ti, ENTRY(dd, i), ctx_rnd_ceil);
                    status |= gr_abs(tj, ENTRY(du, i), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
            }
        }
        gr_heap_clear_vec(dd, n, ctx);
        gr_heap_clear_vec(du, n, ctx);

        gr_mat_init(mm, n, n, ctx);
        gr_mat_init(rm, n, n, ctx);
        gr_mat_init(mn, n, n, ctx);
        gr_mat_init(rn, n, n, ctx);
        gr_mat_init(bound2, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(tmp, MAT_ENTRY(Wu, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_div_ui(MAT_ENTRY(mm, j, i), tmp, 2, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rm, j, i), MAT_ENTRY(mm, j, i), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_div_ui(MAT_ENTRY(mn, i, j), tmp, 2, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rn, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(Wd, mm, mn, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wd, i, i), MAT_ENTRY(Wd, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul_classical(Wu, mm, mn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wu, i, i), MAT_ENTRY(Wu, i, i), 1, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                status |= gr_max(MAT_ENTRY(Wu, i, j), ti, tj, ctx_rnd_ceil);
                status |= gr_abs(MAT_ENTRY(mm, i, j), MAT_ENTRY(mm, i, j), ctx_rnd_ceil);
                status |= gr_abs(MAT_ENTRY(mn, i, j), MAT_ENTRY(mn, i, j), ctx_rnd_ceil);
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(rn, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound2, rm, bound2, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(bound2, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);

            }
        }
        status |= gr_mat_mul(Wu, mm, rn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(bound2, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);

            }
        }

        gr_mat_clear(Wu, ctx);
        gr_mat_clear(Wd, ctx);
        gr_mat_clear(mm, ctx);
        gr_mat_clear(mn, ctx);
        gr_mat_clear(rm, ctx);
        gr_mat_clear(rn, ctx);

        gr_mat_init(Wu, m, n, ctx);
        gr_mat_init(Wd, m, n, ctx);
        gr_mat_init(mm, n, m, ctx);
        gr_mat_init(mn, m, n, ctx);
        gr_mat_init(rm, n, m, ctx);
        gr_mat_init(rn, m, n, ctx);

        status |= gr_mat_mul(Wd, A, V, ctx_rnd_floor);
        status |= gr_mat_mul(Wu, A, V, ctx_rnd_ceil);

        gr_mat_clear(A, ctx);
        gr_mat_clear(V, ctx);

        gr_mat_init(bound3, n, n, ctx);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(tmp, MAT_ENTRY(Wu, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_div_ui(MAT_ENTRY(mm, j, i), tmp, 2, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rm, j, i), MAT_ENTRY(mm, j, i), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_div_ui(MAT_ENTRY(mn, i, j), tmp, 2, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rn, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
            }
        }

        gr_mat_clear(Wd, ctx);
        gr_mat_clear(Wu, ctx);

        gr_mat_init(Wd, n, n, ctx);
        gr_mat_init(Wu, n, n, ctx);

        status |= gr_mat_mul(Wd, mm, mn, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wd, i, i), MAT_ENTRY(Wd, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(Wu, mm, mn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wu, i, i), MAT_ENTRY(Wu, i, i), 1, ctx_rnd_ceil);
            for (j = 0; j < m; j++)
            {
                if (j < n)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(Wu, i, j), ti, tj, ctx_rnd_ceil);
                }
                status |= gr_abs(MAT_ENTRY(mm, i, j), MAT_ENTRY(mm, i, j), ctx_rnd_ceil);
            }
        }

        gr_mat_clear(Wd, ctx);
        gr_mat_init(Wd, m, n, ctx);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(MAT_ENTRY(mn, i, j), MAT_ENTRY(mn, i, j), ctx_rnd_ceil);
                status |= gr_add(MAT_ENTRY(Wd, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(rn, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound3, rm, Wd, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound3, i, j), MAT_ENTRY(bound3, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(Wu, mm, rn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound3, i, j), MAT_ENTRY(bound3, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
            }
        }

        gr_mat_clear(Wu, ctx);
        gr_mat_clear(Wd, ctx);
        gr_mat_clear(mm, ctx);
        gr_mat_clear(mn, ctx);
        gr_mat_clear(rm, ctx);
        gr_mat_clear(rn, ctx);

        gr_mat_init(boundt, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_set(MAT_ENTRY(boundt, j, i), MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                status |= gr_set(ti, MAT_ENTRY(bound2, i, j), ctx_rnd_ceil);
                status |= gr_set(tj, MAT_ENTRY(bound3, i, j), ctx_rnd_ceil);
                status |= gr_add(MAT_ENTRY(bound2, i, j), ti, tj, ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound, bound2, bound, ctx_rnd_ceil);
        status |= gr_mat_mul(bound, boundt, bound, ctx_rnd_ceil);

        gr_mat_clear(bound2, ctx);
        gr_mat_clear(bound3, ctx);
        gr_mat_clear(boundt, ctx);

        status |= gr_zero(norm, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_zero(s, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(tmp, MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                status |= gr_add(s, s, tmp, ctx_rnd_ceil);
            }
            status |= gr_max(norm, norm, s, ctx_rnd_ceil);
        }
        if (status != GR_SUCCESS || gr_lt(norm, one, ctx) != T_TRUE)
            goto fail_clear_R_bound_bla;

        gr_mat_init(absR, n, n, ctx);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j >= i)
                {
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_zero(MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                }
                status |= gr_abs(MAT_ENTRY(absR, i, j), MAT_ENTRY(R, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound, bound, absR, ctx_rnd_ceil);

        gr_mat_clear(absR, ctx);

        for (i = 0; i < n - 1; i++)
        {
            status |= gr_sub(tmp, MAT_ENTRY(R, i, i), MAT_ENTRY(bound, i, i), ctx_rnd_floor);

            /* want  status |= gr_mul_d(ti, tmp, fl->eta, ctx_rnd_floor); */
            status |= gr_set_d(ti, fl->eta, ctx_rnd_ceil);  /* actually exact; rnd does not matter */
            status |= gr_mul(ti, tmp, ti, ctx_rnd_floor);

            for (j = i + 1; j < n; j++)
            {
                status |= gr_abs(tmp, MAT_ENTRY(R, i, j), ctx_rnd_ceil);
                status |= gr_add(tj, tmp, MAT_ENTRY(bound, i, j), ctx_rnd_ceil);

                if (status != GR_SUCCESS || gr_le(tj, ti, ctx) != T_TRUE)
                    goto fail_clear_R_bound_bla;
            }
            status |= gr_add(ti, MAT_ENTRY(R, i, i), MAT_ENTRY(bound, i, i), ctx_rnd_ceil);
            status |= gr_sub(tj, MAT_ENTRY(R, i + 1, i + 1), MAT_ENTRY(bound, i + 1, i + 1), ctx_rnd_floor);
            status |= gr_abs(tmp, MAT_ENTRY(R, i, i + 1), ctx_rnd_floor);
            status |= gr_sub(norm, tmp, MAT_ENTRY(bound, i, i + 1), ctx_rnd_floor);
            status |= gr_div(tmp, norm, ti, ctx_rnd_floor);
            status |= gr_sqr(norm, tmp, ctx_rnd_floor);

            /* want status |= gr_sub_d(s, norm, fl->delta, ctx_rnd_floor); */
            status |= gr_set_d(tmp, fl->delta, ctx_rnd_ceil);  /* actually exact; rnd does not matter */
            status |= gr_sub(s, norm, tmp, ctx_rnd_floor);

            status |= gr_neg(s, s, ctx_rnd_floor);
            status |= gr_sqrt(tmp, s, ctx_rnd_ceil);
            status |= gr_mul(s, tmp, ti, ctx_rnd_ceil);

            if (status != GR_SUCCESS || gr_le(s, tj, ctx) != T_TRUE)
                goto fail_clear_R_bound_bla;
        }

        gr_mat_clear(R, ctx);
        gr_mat_clear(bound, ctx);
        GR_TMP_CLEAR5(s, norm, ti, tj, tmp, ctx);
    }
    else
    {
        slong i, j, k, m;
        gr_mat_t bound2, bound3, boundt, mm, rm, mn, rn, absR;

        if (B->r == 0 || B->r == 1)
        {
            GR_TMP_CLEAR2(zero, one, ctx);
            return 1;
        }

        m = B->c;
        n = B->r;

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(R, n, n, ctx);
        gr_mat_init(V, n, n, ctx);

        GR_TMP_INIT5(s, norm, ti, tj, tmp, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                status |= gr_set_fmpz(MAT_ENTRY(A, j, i), fmpz_mat_entry(B, i, j), ctx);
            }
        }

        /* Todo: implement using dot products */
        for (j = 0; j < n; j++)
        {
            status |= gr_set(MAT_ENTRY(R, j, j), MAT_ENTRY(A, j, j), ctx);
            for (i = 0; i <= j - 1; i++)
            {
                status |= gr_set(MAT_ENTRY(R, i, j), MAT_ENTRY(A, j, i), ctx);
                for (k = 0; k <= i - 1; k++)
                {
                    status |= gr_mul(tmp, MAT_ENTRY(R, k, i), MAT_ENTRY(R, k, j), ctx);
                    status |= gr_sub(MAT_ENTRY(R, i, j), MAT_ENTRY(R, i, j), tmp, ctx);
                }

                if (gr_is_zero(MAT_ENTRY(R, i, i), ctx) != T_TRUE)
                {
                    status |= gr_div(MAT_ENTRY(R, i, j), MAT_ENTRY(R, i, j), MAT_ENTRY(R, i, i), ctx);
                    status |= gr_mul(tmp, MAT_ENTRY(R, i, j), MAT_ENTRY(R, i, j), ctx);
                    status |= gr_sub(MAT_ENTRY(R, j, j), MAT_ENTRY(R, j, j), tmp, ctx);
                }
            }

            if (status != GR_SUCCESS || gr_gt(MAT_ENTRY(R, j, j), zero, ctx) != T_TRUE)
            {
                /* going to take sqrt and then divide by it */
                goto fail_clear_A_R_V;
            }

            status |= gr_sqrt(MAT_ENTRY(R, j, j), MAT_ENTRY(R, j, j), ctx);
        }

        for (j = n - 1; j >= 0; j--)
        {
            status |= gr_inv(MAT_ENTRY(V, j, j), MAT_ENTRY(R, j, j), ctx);
            for (i = j + 1; i < n; i++)
            {
                for (k = j + 1; k < n; k++)
                {
                    status |= gr_mul(norm, MAT_ENTRY(V, k, i), MAT_ENTRY(R, j, k), ctx);
                    status |= gr_add(MAT_ENTRY(V, j, i), MAT_ENTRY(V, j, i), norm, ctx);
                }
                status |= gr_neg(MAT_ENTRY(V, j, i), MAT_ENTRY(V, j, i), ctx);
                status |= gr_mul(MAT_ENTRY(V, j, i), MAT_ENTRY(V, j, i), MAT_ENTRY(V, j, j), ctx);
            }
        }

        gr_mat_init(Wu, n, n, ctx);
        gr_mat_init(Wd, n, n, ctx);
        du = gr_heap_init_vec(n, ctx);
        dd = gr_heap_init_vec(n, ctx);

        status |= gr_mat_mul(Wd, R, V, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(dd, i), MAT_ENTRY(Wd, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(Wu, R, V, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(du, i), MAT_ENTRY(Wu, i, i), 1, ctx_rnd_ceil);
        }
        status |= gr_zero(norm, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_zero(s, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(s, s, tmp, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_abs(ti, ENTRY(dd, i), ctx_rnd_ceil);
                    status |= gr_abs(tj, ENTRY(du, i), ctx_rnd_ceil);
                    status |= gr_max(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(s, s, tmp, ctx_rnd_ceil);
                }
            }
            status |= gr_max(norm, norm, s, ctx_rnd_ceil);
        }

        if (status != GR_SUCCESS || gr_lt(norm, one, ctx) != T_TRUE)
        {
fail_clear_all:
            gr_mat_clear(Wu, ctx);
            gr_mat_clear(Wd, ctx);
            gr_heap_clear_vec(du, n, ctx);
            gr_heap_clear_vec(dd, n, ctx);
            GR_TMP_CLEAR5(s, norm, ti, tj, tmp, ctx);
fail_clear_A_R_V:
            gr_mat_clear(A, ctx);
            gr_mat_clear(R, ctx);
            gr_mat_clear(V, ctx);

            GR_TMP_CLEAR2(zero, one, ctx);
            return 0;
        }

        gr_mat_init(bound, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(dd, i), MAT_ENTRY(Wd, i, i), 2, ctx_rnd_floor);
        }
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(ENTRY(du, i), MAT_ENTRY(Wu, i, i), 2, ctx_rnd_ceil);
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j > i)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
                else if (j < i)
                {
                    status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                    status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_abs(ti, ENTRY(dd, i), ctx_rnd_ceil);
                    status |= gr_abs(tj, ENTRY(du, i), ctx_rnd_ceil);
                    status |= gr_max(MAT_ENTRY(bound, i, j), ti, tj, ctx_rnd_ceil);
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
            }
        }
        gr_heap_clear_vec(dd, n, ctx);
        gr_heap_clear_vec(du, n, ctx);

        gr_mat_init(mm, n, n, ctx);
        gr_mat_init(rm, n, n, ctx);
        gr_mat_init(mn, n, n, ctx);
        gr_mat_init(rn, n, n, ctx);
        gr_mat_init(bound2, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(tmp, MAT_ENTRY(Wu, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_mul_2exp_si(MAT_ENTRY(mm, j, i), tmp, -1, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rm, j, i), MAT_ENTRY(mm, j, i), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_mul_2exp_si(MAT_ENTRY(mn, i, j), tmp, -1, ctx_rnd_ceil);
                status |= gr_sub(MAT_ENTRY(rn, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(Wd, mm, mn, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wd, i, i), MAT_ENTRY(Wd, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(Wu, mm, mn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(Wu, i, i), MAT_ENTRY(Wu, i, i), 1, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(ti, MAT_ENTRY(Wd, i, j), ctx_rnd_ceil);
                status |= gr_abs(tj, MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
                status |= gr_max(MAT_ENTRY(Wu, i, j), ti, tj, ctx_rnd_ceil);
                status |= gr_abs(MAT_ENTRY(mm, i, j), MAT_ENTRY(mm, i, j), ctx_rnd_ceil);
                status |= gr_abs(MAT_ENTRY(mn, i, j), MAT_ENTRY(mn, i, j), ctx_rnd_ceil);
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(mn, i, j), MAT_ENTRY(rn, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound2, rm, bound2, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(bound2, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(Wu, mm, rn, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_add(MAT_ENTRY(bound2, i, j), MAT_ENTRY(bound2, i, j), MAT_ENTRY(Wu, i, j), ctx_rnd_ceil);
            }
        }

        gr_mat_clear(Wu, ctx);
        gr_mat_clear(Wd, ctx);
        gr_mat_clear(mm, ctx);
        gr_mat_clear(mn, ctx);
        gr_mat_clear(rm, ctx);
        gr_mat_clear(rn, ctx);

        gr_mat_init(Wu, m, n, ctx);
        gr_mat_init(Wd, m, n, ctx);
        gr_mat_init(mm, n, m, ctx);
        gr_mat_init(mn, m, n, ctx);
        gr_mat_init(rm, n, m, ctx);
        gr_mat_init(rn, m, n, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_set(MAT_ENTRY(mm, j, i), MAT_ENTRY(V, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(Wd, mm, A, ctx_rnd_floor);
        status |= gr_mat_mul(Wu, mm, A, ctx_rnd_ceil);

        gr_mat_clear(A, ctx);

        gr_mat_init(bound3, n, n, ctx);

        status |= gr_mat_mul(mm, Wd, V, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(mm, i, i), MAT_ENTRY(mm, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(rm, Wd, V, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(rm, i, i), MAT_ENTRY(rm, i, i), 1, ctx_rnd_ceil);
        }

        status |= gr_mat_mul(mn, Wu, V, ctx_rnd_floor);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(mn, i, i), MAT_ENTRY(mn, i, i), 1, ctx_rnd_floor);
        }
        status |= gr_mat_mul(rn, Wu, V, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_sub_ui(MAT_ENTRY(rn, i, i), MAT_ENTRY(rn, i, i), 1, ctx_rnd_ceil);
        }

        gr_mat_clear(Wd, ctx);
        gr_mat_clear(Wu, ctx);
        gr_mat_clear(V, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(ti, MAT_ENTRY(mm, i, j), ctx_rnd_ceil);
                status |= gr_abs(tj, MAT_ENTRY(mn, i, j), ctx_rnd_ceil);
                status |= gr_max(MAT_ENTRY(bound3, i, j), ti, tj, ctx_rnd_ceil);
                status |= gr_abs(tmp, MAT_ENTRY(rm, i, j), ctx_rnd_ceil);
                status |= gr_max(MAT_ENTRY(bound3, i, j), MAT_ENTRY(bound3, i, j), tmp, ctx_rnd_ceil);
                status |= gr_abs(tmp, MAT_ENTRY(rn, i, j), ctx_rnd_ceil);
                status |= gr_max(MAT_ENTRY(bound3, i, j), MAT_ENTRY(bound3, i, j), tmp, ctx_rnd_ceil);
            }
        }

        gr_mat_clear(mm, ctx);
        gr_mat_clear(mn, ctx);
        gr_mat_clear(rm, ctx);
        gr_mat_clear(rn, ctx);

        gr_mat_init(boundt, n, n, ctx);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_set(MAT_ENTRY(boundt, j, i), MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                status |= gr_set(ti, MAT_ENTRY(bound2, i, j), ctx_rnd_ceil);
                status |= gr_set(tj, MAT_ENTRY(bound3, i, j), ctx_rnd_ceil);
                status |= gr_add(MAT_ENTRY(bound2, i, j), ti, tj, ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound, bound2, bound, ctx_rnd_ceil);
        status |= gr_mat_mul(bound, boundt, bound, ctx_rnd_ceil);

        gr_mat_clear(bound2, ctx);
        gr_mat_clear(bound3, ctx);
        gr_mat_clear(boundt, ctx);

        status |= gr_zero(norm, ctx_rnd_ceil);
        for (i = 0; i < n; i++)
        {
            status |= gr_zero(s, ctx_rnd_ceil);
            for (j = 0; j < n; j++)
            {
                status |= gr_abs(tmp, MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                status |= gr_add(s, s, tmp, ctx_rnd_ceil);
            }
            status |= gr_max(norm, norm, s, ctx_rnd_ceil);
        }
        if (status != GR_SUCCESS || gr_lt(norm, one, ctx) != T_TRUE)
            goto fail_clear_R_bound_bla;


        gr_mat_init(absR, n, n, ctx);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j >= i)
                {
                    status |= gr_mul(ti, norm, norm, ctx_rnd_ceil);
                    status |= gr_sub(tj, one, norm, ctx_rnd_ceil);
                    status |= gr_div(tmp, ti, tj, ctx_rnd_ceil);
                    status |= gr_add(MAT_ENTRY(bound, i, j), MAT_ENTRY(bound, i, j), tmp, ctx_rnd_ceil);
                }
                else
                {
                    status |= gr_zero(MAT_ENTRY(bound, i, j), ctx_rnd_ceil);
                }
                status |= gr_abs(MAT_ENTRY(absR, i, j), MAT_ENTRY(R, i, j), ctx_rnd_ceil);
            }
        }
        status |= gr_mat_mul(bound, bound, absR, ctx_rnd_ceil);

        gr_mat_clear(absR, ctx);

        for (i = 0; i < n - 1; i++)
        {
            status |= gr_sub(tmp, MAT_ENTRY(R, i, i), MAT_ENTRY(bound, i, i), ctx_rnd_floor);

            /* want  status |= gr_mul_d(ti, tmp, fl->eta, ctx_rnd_floor); */
            status |= gr_set_d(ti, fl->eta, ctx_rnd_ceil);  /* actually exact; rnd does not matter */
            status |= gr_mul(ti, tmp, ti, ctx_rnd_floor);

            for (j = i + 1; j < n; j++)
            {
                status |= gr_abs(tmp, MAT_ENTRY(R, i, j), ctx_rnd_ceil);
                status |= gr_add(tj, tmp, MAT_ENTRY(bound, i, j), ctx_rnd_ceil);

                if (status != GR_SUCCESS || gr_le(tj, ti, ctx) != T_TRUE)
                    goto fail_clear_R_bound_bla;
            }
            status |= gr_add(ti, MAT_ENTRY(R, i, i), MAT_ENTRY(bound, i, i), ctx_rnd_ceil);
            status |= gr_sub(tj, MAT_ENTRY(R, i + 1, i + 1), MAT_ENTRY(bound, i + 1, i + 1), ctx_rnd_floor);
            status |= gr_abs(tmp, MAT_ENTRY(R, i, i + 1), ctx_rnd_floor);
            status |= gr_sub(norm, tmp, MAT_ENTRY(bound, i, i + 1), ctx_rnd_floor);
            status |= gr_div(tmp, norm, ti, ctx_rnd_floor);
            status |= gr_sqr(norm, tmp, ctx_rnd_floor);

            /* want status |= gr_sub_d(s, norm, fl->delta, ctx_rnd_floor); */
            status |= gr_set_d(tmp, fl->delta, ctx_rnd_ceil);  /* actually exact; rnd does not matter */
            status |= gr_sub(s, norm, tmp, ctx_rnd_floor);

            status |= gr_neg(s, s, ctx_rnd_floor);
            status |= gr_sqrt(tmp, s, ctx_rnd_ceil);
            status |= gr_mul(s, tmp, ti, ctx_rnd_ceil);

            if (status != GR_SUCCESS || gr_le(s, tj, ctx) != T_TRUE)
            {
fail_clear_R_bound_bla:
                gr_mat_clear(R, ctx);
                gr_mat_clear(bound, ctx);
                GR_TMP_CLEAR5(s, norm, ti, tj, tmp, ctx);

                GR_TMP_CLEAR2(zero, one, ctx);
                return 0;
            }
        }

        gr_mat_clear(R, ctx);
        gr_mat_clear(bound, ctx);
        GR_TMP_CLEAR5(s, norm, ti, tj, tmp, ctx);
    }

    FLINT_ASSERT((fl->rt == Z_BASIS
                        ? fmpz_mat_is_reduced(B, fl->delta, fl->eta)
                        : fmpz_mat_is_reduced_gram(B, fl->delta, fl->eta)));

    GR_TMP_CLEAR2(zero, one, ctx);

    return 1;
}

