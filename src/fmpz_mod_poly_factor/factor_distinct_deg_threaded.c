/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <math.h>
#include <pthread.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "fmpz_mod_poly.h"
#include "thread_support.h"

void
_fmpz_mod_poly_precompute_matrix_worker (void * arg_ptr)
{
    fmpz_mod_poly_matrix_precompute_arg_t arg =
                           *((fmpz_mod_poly_matrix_precompute_arg_t *) arg_ptr);
    /* Set rows of A to powers of poly1 */
    slong i, n, m;
    fmpz_mod_poly_struct * poly1 = arg.poly1;
    fmpz_mod_poly_struct * poly2 = arg.poly2;
    fmpz_mod_poly_struct * poly2inv = arg.poly2inv;
    fmpz_mat_struct * A = arg.A;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    n = poly2->length - 1;

    m = n_sqrt(n) + 1;

    fmpz_one(A->rows[0] + 0);
    _fmpz_vec_set(A->rows[1], poly1->coeffs, n);

    for (i = 2; i < m; i++)
        _fmpz_mod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n,
                                     poly1->coeffs, n, poly2->coeffs,
                                     n + 1, poly2inv->coeffs, n + 1, p);
}

void
_fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr)
{
    fmpz_mod_poly_compose_mod_precomp_preinv_arg_t arg=
                   *((fmpz_mod_poly_compose_mod_precomp_preinv_arg_t*) arg_ptr);
    fmpz_mat_t B, C;
    fmpz * t, * h;
    slong i, j, n, m;
    fmpz_mod_poly_struct * res = arg.res;
    fmpz_mod_poly_struct * poly1 = arg.poly1;
    fmpz_mod_poly_struct * poly3 = arg.poly3;
    fmpz_mod_poly_struct * poly3inv = arg.poly3inv;
    fmpz_mat_struct * A = arg.A;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    if (poly3->length == 1)
        return;

    if (poly1->length == 1)
    {
        fmpz_set(res->coeffs, poly1->coeffs);
        return;
    }

    if (poly3->length == 2)
    {
        _fmpz_mod_poly_evaluate_fmpz(res->coeffs, poly1->coeffs,
                                poly1->length, A->rows[1] + 0, p);
        return;
    }

    n = poly3->length - 1;
    m = n_sqrt(n) + 1;

    fmpz_mat_init(B, m, m);
    fmpz_mat_init(C, m, n);

    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < poly1->length / m; i++)
        _fmpz_vec_set(B->rows[i], poly1->coeffs + i*m, m);

    _fmpz_vec_set(B->rows[i], poly1->coeffs + i*m, poly1->length % m);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);
    }

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res->coeffs, C->rows[m - 1], n);
    _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n,
                                 poly3->coeffs, poly3->length,
                                 poly3inv->coeffs, poly3inv->length, p);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_mulmod_preinv(t, res->coeffs, n, h, n,
                                     poly3->coeffs, poly3->length,
                                     poly3inv->coeffs, poly3->length, p);
        _fmpz_mod_poly_add(res->coeffs, t, n, C->rows[i], n, p);
    }

    _fmpz_vec_clear(h, n);
    _fmpz_vec_clear(t, n);

    fmpz_mat_clear(B);
    fmpz_mat_clear(C);

    return;
}

void
_fmpz_mod_poly_interval_poly_worker(void * arg_ptr)
{
    fmpz_mod_poly_interval_poly_arg_t arg =
                               *((fmpz_mod_poly_interval_poly_arg_t *) arg_ptr);
    slong k, m = arg.m;
    fmpz_mod_poly_struct * H = arg.H;
    fmpz_mod_poly_struct * res = arg.res;
    fmpz_mod_poly_struct * v = arg.v;
    fmpz_mod_poly_struct * vinv = arg.vinv;
    fmpz_mod_poly_struct * baby = arg.baby;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz * tmp = arg.tmp;
    fmpz_t invV;

    fmpz_init(invV);

    fmpz_invmod(invV, v->coeffs + v->length - 1, p);

    fmpz_set_ui(res->coeffs + 0, 1);

    for (k = m - 1; k >= 0; k--)
    {
        _fmpz_vec_zero(tmp, v->length - 1);

        if (baby[k].length < v->length)
          _fmpz_vec_set(tmp, baby[k].coeffs, baby[k].length);
        else
          _fmpz_mod_poly_rem(tmp, baby[k].coeffs, baby[k].length,
                             v->coeffs , v->length, invV, p);

        _fmpz_mod_poly_sub(tmp, H->coeffs, H->length, tmp,
                           v->length - 1, p);

        _fmpz_mod_poly_mulmod_preinv(res->coeffs, tmp, v->length - 1,
                                     res->coeffs, v->length - 1,
                                     v->coeffs, v->length,
                                     vinv->coeffs, vinv->length, p);
    }

    fmpz_clear(invV);

    return;
}

/* the degrees are written as exponents of the corresponding factors */
void fmpz_mod_poly_factor_distinct_deg_threaded_with_frob(
    fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t poly,
    const fmpz_mod_poly_t polyinv,
    const fmpz_mod_poly_t frob,  /* x^p mod poly */
    const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t f, g, v, vinv, tmp, II;
    fmpz_mod_poly_struct * h, * H, * I, * scratch;
    slong i, j, k, l, m, n, d, c1 = 1, c2;
    fmpz_mat_struct * HH;
    double beta;
    thread_pool_handle * threads;
    slong num_threads;
    fmpz_mod_poly_matrix_precompute_arg_t * args1;
    fmpz_mod_poly_compose_mod_precomp_preinv_arg_t * args2;
    fmpz_mod_poly_interval_poly_arg_t * args3;

    FLINT_ASSERT(fmpz_mod_poly_is_monic(poly, ctx));

    n = fmpz_mod_poly_degree(poly, ctx);
    if (n == 1)
    {
        fmpz_mod_poly_factor_fit_length(res, 1, ctx);
        fmpz_mod_poly_set(res->poly + 0, poly, ctx);
        res->exp[0] = 1;
        res->num = 1;
        return;
    }

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */
    fmpz_mod_poly_init(v, ctx);
    fmpz_mod_poly_init(vinv, ctx);
    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(tmp, ctx);
    fmpz_mod_poly_init(II, ctx);

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    h = FLINT_ARRAY_ALLOC(2*m + l + num_threads + 2, fmpz_mod_poly_struct);

    for (i = 0; i < 2*m + l + 2 + num_threads; i++)
       fmpz_mod_poly_init(h + i, ctx);

    H = h + (l + 1);
    I = H + m;
    scratch = I + m;

    HH    = FLINT_ARRAY_ALLOC(num_threads + 2, fmpz_mat_struct);
    args1 = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mod_poly_matrix_precompute_arg_t);
    args2 = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mod_poly_compose_mod_precomp_preinv_arg_t);
    args3 = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mod_poly_interval_poly_arg_t);

    fmpz_mod_poly_set(v, poly, ctx);
    fmpz_mod_poly_set(vinv, polyinv, ctx);

    /* compute baby steps: h[i]=x^{p^i}mod v */
    fmpz_mod_poly_set_coeff_ui(h + 0, 1, 1, ctx);
    fmpz_mod_poly_set(h + 1, frob, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_powmod_x_fmpz_preinv(tmp, p, v, vinv, ctx);
    FLINT_ASSERT(fmpz_mod_poly_equal(tmp, h + 1, ctx));
#endif

    if (fmpz_sizeinbase(p, 2) > ((n_sqrt(v->length - 1) + 1) * 3) / 4)
    {
        for (i = 1; i < FLINT_BIT_COUNT(l); i++)
            fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(h + 1 +
                                                             (1 << (i - 1)),
                                                             h + 1,
                                                             1 << (i - 1),
                                                             1 << (i - 1),
                                                             h + (1 << (i - 1)),
                                                             v, vinv, ctx,
                                                             threads, num_threads);

        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(h + 1 +
                                                         (1 << (i - 1)),
                                                         h + 1,
                                                         1 << (i - 1),
                                                         l - (1 << (i - 1)),
                                                         h + (1 << (i - 1)),
                                                         v, vinv, ctx,
                                                         threads, num_threads);
    } else
    {
        for (i = 2; i < l + 1; i++)
        {
            fmpz_mod_poly_init(h + i, ctx);
            fmpz_mod_poly_powmod_fmpz_binexp_preinv(h + i, h + i - 1, p,
                                                                 v, vinv, ctx);
        }
    }

    /* compute coarse distinct-degree factorisation */
    res->num = 0;
    fmpz_mod_poly_set(H + 0, h + l, ctx);
    fmpz_mat_init(HH + 0, n_sqrt(v->length - 1) + 1, v->length - 1);

    fmpz_mod_poly_precompute_matrix(HH + 0, H + 0, v, vinv, ctx);

    for (d = 1, j = 0; j < m /(num_threads + 1) + 1; j++)
    {
        if (j == 0)
        {
            for (i = 0; i < num_threads + 1; i++)
            {
                if (i > 0 && I[i - 1].length > 1)
                {
                    _fmpz_mod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                          HH + 0, v, ctx);

                    fmpz_mat_clear(HH + 0);
                    fmpz_mat_init_set(HH + 0, HH + num_threads + 1);

                    fmpz_mat_clear(HH + num_threads + 1);
                    
                    fmpz_mod_poly_rem(tmp, H + i - 1, v, ctx);

                    fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H + i,
                                                    tmp, HH + 0, v, vinv, ctx);
                } else if (i > 0)
                    fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H + i,
                                              H + i - 1, HH + 0, v, vinv, ctx);

                /* compute interval polynomials */
                fmpz_mod_poly_set_coeff_ui(I + i, 0, 1, ctx);

                for (k = l - 1; k >= 0 && 2*d <= v->length - 1; k--, d++)
                {
                    fmpz_mod_poly_rem(tmp, h + k, v, ctx);
                    fmpz_mod_poly_sub(tmp, H + i, tmp, ctx);
                    fmpz_mod_poly_mulmod_preinv(I + i, tmp, I + i, v, vinv, ctx);
                }

                /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
                /* F_j is stored on the place of I_j */
                fmpz_mod_poly_gcd(I + i, v, I + i, ctx);

                if (I[i].length > 1)
                {
                    fmpz_mod_poly_remove(v, I + i, ctx);
                    fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
                    fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);
                }

                if (v->length - 1 < 2*d)
                    break;
            }

            if (v->length - 1 < 2*d)
                break;
        } else if (j == 1 && num_threads + 1 < m)
        {
            if (I[num_threads].length > 1)
            {
                _fmpz_mod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                               HH + 0, v, ctx);

                fmpz_mat_clear(HH + 0);
                fmpz_mat_init_set(HH + 0, HH + num_threads + 1);
                fmpz_mat_clear(HH + num_threads + 1);
            }

            /* we make thread 0 the master thread, but don't use it in this loop */
            for (c1 = 1, i = 1; i < num_threads + 1 &&
                                    i + num_threads + 1 < m; i++, c1++)
            {
                fmpz_mat_init(HH + i, n_sqrt(v->length - 1) + 1, v->length - 1);

                fmpz_mod_poly_rem(scratch + i, H + i, v, ctx);

                if (scratch[i].length < v->length - 1)
                {
                    fmpz_mod_poly_fit_length(scratch + i, v->length - 1, ctx);
                    _fmpz_vec_zero(scratch[i].coeffs + scratch[i].length,
                                   v->length - 1 - scratch[i].length);
                    _fmpz_mod_poly_set_length(scratch + i, v->length - 1);
                }

                args1[i].A        = HH + i;
                args1[i].poly1    = scratch + i;
                args1[i].poly2    = v;
                args1[i].poly2inv = vinv;
                args1[i].ctx      = ctx;

                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                            _fmpz_mod_poly_precompute_matrix_worker, &args1[i]);
            }
            
            for (i = 1; i < c1; i++)
                thread_pool_wait(global_thread_pool, threads[i - 1]);

            fmpz_mod_poly_rem(tmp, H + num_threads, v, ctx);
            
            for (i = 0; i < c1; i++)
            {
                fmpz_mod_poly_fit_length(H + num_threads + 1 + i, v->length - 1, ctx);
                _fmpz_mod_poly_set_length(H + num_threads + 1 + i, v->length - 1);
                
                args2[i].A        = HH + i;
                args2[i].res      = H + num_threads + i + 1;
                args2[i].poly1    = tmp;
                args2[i].poly3    = v;
                args2[i].poly3inv = vinv;
                args2[i].ctx      = ctx;
            }

            for (i = 1; i < c1; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
        _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args2[i]);
            }
            
            _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args2[0]);
            
            _fmpz_mod_poly_normalise(H + num_threads + 1);
            
            for (i = 1; i < c1; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _fmpz_mod_poly_normalise(H + num_threads + 1 + i);
            }

            for (i = 0; i < c1; i++)
            {
                fmpz_mod_poly_fit_length(I + num_threads + 1 + i, v->length - 1, ctx);
                _fmpz_mod_poly_set_length(I + num_threads + 1 + i, v->length - 1);
                _fmpz_vec_zero(I[num_threads + 1 + i].coeffs, v->length - 1);
                
                args3[i].baby = h;
                args3[i].H    = H + num_threads + 1 + i;
                args3[i].m    = l;
                args3[i].res  = I + num_threads + 1 + i;
                args3[i].v    = v;
                args3[i].vinv = vinv;
                args3[i].ctx  = ctx;
                args3[i].tmp  = _fmpz_vec_init(v->length - 1);
            }

            for (i = 1; i < c1; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                               _fmpz_mod_poly_interval_poly_worker, &args3[i]);
            }

            _fmpz_mod_poly_interval_poly_worker(&args3[0]);
            
            _fmpz_mod_poly_normalise(I + num_threads + 1);

            for (i = 1; i < c1; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _fmpz_mod_poly_normalise(I + num_threads + i + 1);
            }

            for (i = 0; i < c1; i++)
            {
               _fmpz_vec_clear(args3[i].tmp, v->length - 1);
            }

            fmpz_mod_poly_set_ui(II, 1, ctx);

            for (i = 0; i < c1; i++)
                fmpz_mod_poly_mulmod_preinv(II, II,  I + num_threads + 1 + i,
                                                                 v, vinv, ctx);

            fmpz_mod_poly_gcd(II, v, II, ctx);

            if (II->length > 1)
            {
                fmpz_mod_poly_remove(v, II, ctx);
                fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
                fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);
                
                for (i = 0; i < c1; i++)
                {
                    fmpz_mod_poly_gcd(I + num_threads + 1 + i,
                                             I + num_threads + 1 + i, II, ctx);
                    if (I[num_threads + 1 + i].length > 1)
                        fmpz_mod_poly_remove(II, I + num_threads + 1 + i, ctx);
                }
            } else
            {
                for (i = 0; i < c1; i++)
                    fmpz_mod_poly_set_ui(I + num_threads + 1 + i, 1, ctx);
            }

            d = d + c1*l;

            if (v->length-1 < 2*d)
                break;
        } else if (j*(num_threads + 1) < m)
        {
            for (c2 = 0, i = 0; i < num_threads + 1 &&
                                        j*(num_threads + 1) + i < m; i++, c2++)
            {
                if (HH[i].c > v -> length - 1)
                {
                    _fmpz_mod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                               HH + i, v, ctx);
                    
                    fmpz_mat_clear(HH + i);
                    fmpz_mat_init_set(HH + i, HH + num_threads + 1);
                    fmpz_mat_clear(HH + num_threads + 1);
                }
            }

            fmpz_mod_poly_rem(tmp, H + j*(num_threads + 1) - 1, v, ctx);
            
            for (i = 0; i < c2; i++)
            {
                fmpz_mod_poly_fit_length(H + j*(num_threads + 1) + i,
                                                           v->length - 1, ctx);
                _fmpz_mod_poly_set_length(H + j*(num_threads + 1) + i,
                                                                v->length - 1);

                args2[i].A        = HH + i;
                args2[i].res      = H + j*(num_threads + 1) + i;
                args2[i].poly1    = tmp;
                args2[i].poly3    = v;
                args2[i].poly3inv = vinv;
                args2[i].ctx      = ctx;
            }
            
            for (i = 1; i < c2; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
        _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args2[i]);
            }

            _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args2[0]);
            _fmpz_mod_poly_normalise(H + j*(num_threads + 1));

            for (i = 1; i < c2; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _fmpz_mod_poly_normalise(H + j*(num_threads + 1) + i);
            }

            for (i = 0; i < c2; i++)
            {
                fmpz_mod_poly_fit_length(I + j*(num_threads + 1) + i,
                                                           v->length - 1, ctx);
                _fmpz_mod_poly_set_length(I + j*(num_threads + 1) + i,
                                                                v->length - 1);
                _fmpz_vec_zero(I[j*(num_threads + 1) + i].coeffs,
                                                                v->length - 1);
                
                args3[i].baby = h;
                args3[i].H    = H + j*(num_threads + 1) + i;
                args3[i].m    = l;
                args3[i].res  = I + j*(num_threads + 1) + i;
                args3[i].v    = v;
                args3[i].vinv = vinv;
                args3[i].ctx  = ctx;
                args3[i].tmp  = _fmpz_vec_init(v->length - 1);
            }

            for (i = 1; i < c2; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                               _fmpz_mod_poly_interval_poly_worker, &args3[i]);
            }

            _fmpz_mod_poly_interval_poly_worker(&args3[0]);
            _fmpz_mod_poly_normalise(I + j*(num_threads + 1));

            for (i = 1; i < c2; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _fmpz_mod_poly_normalise(I + j*(num_threads + 1) + i);
            }

            for (i = 0; i < c2; i++)
            {
               _fmpz_vec_clear(args3[i].tmp, v->length - 1);
            }

            fmpz_mod_poly_set_ui(II, 1, ctx);

            for (i = 0; i < c2; i++)
                fmpz_mod_poly_mulmod_preinv(II, II,
                                    I + j*(num_threads + 1) + i, v, vinv, ctx);

            fmpz_mod_poly_gcd(II, v, II, ctx);

            if (II->length > 1)
            {
                fmpz_mod_poly_remove(v, II, ctx);

                fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
                fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);

                for (i = 0; i < c2; i++)
                {
                    fmpz_mod_poly_gcd(I + j*(num_threads + 1) + i,
                                      I + j*(num_threads + 1) + i, II, ctx);

                    if (I[j*(num_threads + 1) + i].length > 1)
                        fmpz_mod_poly_remove(II, I + j*(num_threads + 1) + i, ctx);
                }
            } else
            {
                for (i = 0; i < c2; i++)
                    fmpz_mod_poly_set_ui(I + j*(num_threads + 1) + i, 1, ctx);
            }

            d = d + c2*l;

            if (v->length - 1 < 2*d)
                break;
        }
    }

    flint_give_back_threads(threads, num_threads);

    if (v->length > 1)
    {
        fmpz_mod_poly_factor_fit_length(res, res->num + 1, ctx);
        res->exp[res->num] = v->length - 1;
        fmpz_mod_poly_swap(res->poly + res->num, v, ctx);
        res->num++;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j].length - 1 > (j + 1)*l || j == 0)
        {
            fmpz_mod_poly_set(g, I + j, ctx);

            for (i = l - 1; i >= 0 && g->length > 1; i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                fmpz_mod_poly_sub(tmp, H + j, h + i, ctx);
                fmpz_mod_poly_gcd(f, g, tmp, ctx);

                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    fmpz_mod_poly_divrem(g, tmp, g, f, ctx);

                    FLINT_ASSERT(fmpz_mod_poly_is_monic(f, ctx));
                    fmpz_mod_poly_factor_fit_length(res, res->num + 1, ctx);
                    res->exp[res->num] = l * (j + 1) - i;
                    fmpz_mod_poly_swap(res->poly + res->num, f, ctx);
                    res->num++;
                }
            }
        }
        else if (I[j].length > 1)
        {
            FLINT_ASSERT(fmpz_mod_poly_is_monic(I + j, ctx));
            fmpz_mod_poly_factor_fit_length(res, res->num + 1, ctx);
            res->exp[res->num] = I[j].length - 1;
            fmpz_mod_poly_swap(res->poly + res->num, I + j, ctx);
            res->num++;
        }
    }

    /* cleanup */
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_clear(vinv, ctx);
    fmpz_mod_poly_clear(tmp, ctx);
    fmpz_mod_poly_clear(II, ctx);
    
    for (i = 0; i < 2*m + l + num_threads + 2; i++)
        fmpz_mod_poly_clear(h + i, ctx);

    for (i = 0; i < c1; i++)
        fmpz_mat_clear(HH + i);

    flint_free(h);
    flint_free(HH);
    flint_free(args1);
    flint_free(args2);
    flint_free(args3);
}

void fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t poly, slong * const *degs,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t v, vinv, xp;

    fmpz_mod_poly_init(v, ctx);
    fmpz_mod_poly_init(vinv, ctx);
    fmpz_mod_poly_init(xp, ctx);

    fmpz_mod_poly_make_monic(v, poly, ctx);

    fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
    fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);
    fmpz_mod_poly_powmod_x_fmpz_preinv(xp, fmpz_mod_ctx_modulus(ctx), v, vinv, ctx);

    fmpz_mod_poly_factor_distinct_deg_threaded_with_frob(res, v, vinv, xp, ctx);

    /* satisfy the requirements of the interface */
    for (i = 0; i < res->num; i++)
    {
        (*degs)[i] = res->exp[i];
        res->exp[i] = 1;
    }

    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_clear(vinv, ctx);
    fmpz_mod_poly_clear(xp, ctx);
}
