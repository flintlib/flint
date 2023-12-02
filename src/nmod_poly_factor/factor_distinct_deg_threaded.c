/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020, 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
_nmod_poly_precompute_matrix_worker(void * arg_ptr)
{
    nmod_poly_matrix_precompute_arg_t arg =
                           *((nmod_poly_matrix_precompute_arg_t *) arg_ptr);

    /* Set rows of A to powers of poly1 */
    slong i, n, m;
    nmod_poly_struct * poly1 = arg.poly1;
    nmod_poly_struct * poly2 = arg.poly2;
    nmod_poly_struct * poly2inv = arg.poly2inv;
    nmod_mat_struct * A = arg.A;
    nmod_t mod = poly2->mod;

    n = poly2->length - 1;
    m = n_sqrt(n) + 1;

    for (i = 1; i < n; i++)
        A->rows[0][i] = 0;
    A->rows[0][0] = 1;

    _nmod_vec_set(A->rows[1], poly1->coeffs, n);

    for (i = 2; i < m; i++)
        _nmod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n,
                                 poly1->coeffs, n, poly2->coeffs, n + 1,
                                 poly2inv->coeffs, n + 1, mod);
}

void
_nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr)
{
    nmod_poly_compose_mod_precomp_preinv_arg_t arg =
                   *((nmod_poly_compose_mod_precomp_preinv_arg_t*) arg_ptr);
    nmod_mat_t B, C;
    mp_ptr t, h;
    slong i, n, m;
    nmod_poly_struct * res = arg.res;
    nmod_poly_struct * poly1 = arg.poly1;
    nmod_poly_struct * poly3 = arg.poly3;
    nmod_poly_struct * poly3inv = arg.poly3inv;
    nmod_mat_struct * A = arg.A;
    nmod_t mod = poly3->mod;

    if (poly3->length == 1)
        return;

    if (poly1->length == 1)
    {
        res->coeffs[0] = poly1->coeffs[0];
        return;
    }

    if (poly3->length == 2)
    {
        res->coeffs[0] = _nmod_poly_evaluate_nmod(poly1->coeffs,
                                             poly1->length, A->rows[1][0], mod);
        return;
    }

    n = poly3->length - 1;
    m = n_sqrt(n) + 1;

    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < poly1->length/m; i++)
        _nmod_vec_set(B->rows[i], poly1->coeffs + i * m, m);

    _nmod_vec_set(B->rows[i], poly1->coeffs + i * m, poly1->length % m);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res->coeffs, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n,
        poly3->coeffs, poly3->length, poly3inv->coeffs, poly3inv->length, mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, res->coeffs, n, h, n, poly3->coeffs,
                          poly3->length, poly3inv->coeffs, poly3->length, mod);
        _nmod_poly_add(res->coeffs, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
_nmod_poly_interval_poly_worker(void * arg_ptr)
{
    nmod_poly_interval_poly_arg_t arg =
                               *((nmod_poly_interval_poly_arg_t *) arg_ptr);
    slong k, m = arg.m;
    nmod_poly_struct * H = arg.H;
    nmod_poly_struct * res = arg.res;
    nmod_poly_struct * v = arg.v;
    nmod_poly_struct * vinv = arg.vinv;
    nmod_poly_struct * baby = arg.baby;
    nmod_t mod = v->mod;
    mp_ptr tmp = arg.tmp;

    res->coeffs[0] = 1;

    for (k = m - 1; k >= 0; k--)
    {
        flint_mpn_zero(tmp, v->length - 1);

        if (baby[k].length < v->length)
            _nmod_vec_set(tmp, baby[k].coeffs, baby[k].length);
        else
            _nmod_poly_rem(tmp, baby[k].coeffs, baby[k].length,
                                        v->coeffs, v->length, mod);

        _nmod_poly_sub(tmp, H->coeffs, H->length, tmp, v->length - 1, mod);

        _nmod_poly_mulmod_preinv(res->coeffs, tmp, v->length - 1,
                                 res->coeffs, v->length - 1,
                                 v->coeffs, v->length,
                                 vinv->coeffs, vinv->length, mod);
    }
}

void nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res,
                                   const nmod_poly_t poly, slong * const * degs)
{
    nmod_poly_t f, g, v, vinv, tmp, II;
    nmod_poly_struct * h, * H, * I, * scratch;
    slong i, j, k, l, m, n, index, d, c1 = 1, c2;
    nmod_mat_struct * HH;
    double beta;
    thread_pool_handle * threads;
    slong num_threads;
    nmod_poly_matrix_precompute_arg_t * args1;
    nmod_poly_compose_mod_precomp_preinv_arg_t * args2;
    nmod_poly_interval_poly_arg_t * args3;

    n = nmod_poly_degree(poly);
    nmod_poly_init_mod(v, poly->mod);

    nmod_poly_make_monic(v, poly);

    if (n == 1)
    {
        nmod_poly_factor_insert(res, v, 1);
        (*degs)[0] = 1;

        nmod_poly_clear(v);

        return;
    }

    beta = 0.5 * (1. - log(2)/log(n));
    l = ceil(pow(n, beta));
    m = ceil(0.5*n/l);

    /* initialization */
    nmod_poly_init_mod(f, poly->mod);
    nmod_poly_init_mod(g, poly->mod);
    nmod_poly_init_mod(vinv, poly->mod);
    nmod_poly_init_mod(tmp, poly->mod);
    nmod_poly_init_mod(II, poly->mod);

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    h = flint_malloc((2 * m + l + num_threads + 2) * sizeof(nmod_poly_struct));

    for (i = 0; i < 2*m + l + 2 + num_threads; i++)
        nmod_poly_init_mod(h + i, poly->mod);

    H = h + (l + 1);
    I = H + m;
    scratch = I + m;

    HH      = (nmod_mat_struct *)
	               flint_malloc(sizeof(nmod_mat_struct)*(num_threads + 2));
    args1   = (nmod_poly_matrix_precompute_arg_t *)
	               flint_malloc((num_threads + 1)*
                           sizeof(nmod_poly_matrix_precompute_arg_t));
    args2   = (nmod_poly_compose_mod_precomp_preinv_arg_t *)
	               flint_malloc((num_threads + 1)*
                           sizeof(nmod_poly_compose_mod_precomp_preinv_arg_t));
    args3   = (nmod_poly_interval_poly_arg_t *)
	               flint_malloc((num_threads + 1)*
                           sizeof(nmod_poly_interval_poly_arg_t));

    nmod_poly_reverse(vinv, v, v->length);
    nmod_poly_inv_series(vinv, vinv, v->length);

    /* compute baby steps: h[i] = x^{p^i} mod v */
    nmod_poly_set_coeff_ui(h + 0, 1, 1);
    nmod_poly_powmod_x_ui_preinv(h + 1, poly->mod.n, v, vinv);

    if (FLINT_BIT_COUNT(poly->mod.n) > ((n_sqrt(v->length - 1) + 1)*3)/4)
    {
        for (i = 1; i < FLINT_BIT_COUNT(l); i++)
            nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(h + 1 +
                                   (1 << (i - 1)), h + 1, 1 << (i - 1),
                                   1 << (i - 1), h + (1 << (i - 1)),
                                   v, vinv, threads, num_threads);

        nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(h + 1 +
                                   (1 << (i - 1)), h + 1, 1 << (i - 1),
                                   l - (1 << (i - 1)), h + (1 << (i - 1)),
                                   v, vinv, threads, num_threads);
    } else
    {
        for (i = 2; i < l + 1; i++)
        {
            nmod_poly_init_mod(h + i, poly->mod);

            nmod_poly_powmod_ui_binexp_preinv(h + i, h + i - 1, poly->mod.n,
                                                                      v, vinv);
        }
    }

    /* compute coarse distinct-degree factorisation */
    index = 0;

    nmod_poly_set(H + 0, h + l);
    nmod_mat_init(HH + 0, n_sqrt(v->length - 1) + 1,
                                                   v->length - 1, poly->mod.n);

    nmod_poly_precompute_matrix(HH + 0, H + 0, v, vinv);

    for (d = 1, j = 0; j < m/(num_threads + 1) + 1; j++)
    {
        if (j == 0)
        {
            for (i = 0; i < num_threads + 1; i++)
            {
                if (i > 0 && I[i - 1].length > 1)
                {
                    _nmod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                                    HH + 0, v);

                    nmod_mat_clear(HH + 0);
                    nmod_mat_init_set(HH + 0, HH + num_threads + 1);

                    nmod_mat_clear(HH + num_threads + 1);

                    nmod_poly_rem(tmp, H + i - 1, v);

                    nmod_poly_compose_mod_brent_kung_precomp_preinv(H + i,
                                                         tmp, HH + 0, v, vinv);
                } else if (i > 0)
                    nmod_poly_compose_mod_brent_kung_precomp_preinv(H + i,
                                                   H + i - 1, HH + 0, v, vinv);

                /* compute interval polynomials */
                nmod_poly_one(I + i);

                for (k = l - 1; k >= 0 && 2*d <= v->length - 1; k--, d++)
                {
                    nmod_poly_rem(tmp, h + k, v);
                    nmod_poly_sub(tmp, H + i, tmp);
                    nmod_poly_mulmod_preinv(I + i, tmp, I + i, v, vinv);
                }

                /*
                    compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]}
                    F_j is stored on the place of I_j
                */
                nmod_poly_gcd(I + i, v, I + i);

                if (I[i].length > 1)
                {
                    nmod_poly_remove(v, I + i);
                    nmod_poly_reverse(vinv, v, v->length);
                    nmod_poly_inv_series_newton(vinv, vinv, v->length);
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
                _nmod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                                    HH + 0, v);

                nmod_mat_clear(HH + 0);
                nmod_mat_init_set(HH + 0, HH + num_threads + 1);
                nmod_mat_clear(HH + num_threads + 1);
            }

            /* we make thread 0 the master thread, but don't use it in this loop */
            for (c1 = 1, i = 1; i < num_threads + 1 &&
                                    i + num_threads + 1 < m; i++, c1++)
            {
                nmod_mat_init(HH + i, n_sqrt(v->length - 1) + 1, v->length - 1,
                              poly->mod.n);

                nmod_poly_rem(scratch + i, H + i, v);

                if (scratch[i].length < v->length - 1)
                {
                    nmod_poly_fit_length(scratch + i, v->length - 1);
                    flint_mpn_zero(scratch[i].coeffs + scratch[i].length,
                                   v->length - 1 - scratch[i].length);
                    _nmod_poly_set_length(scratch + i, v->length - 1);
                }

                args1[i].A        = HH + i;
                args1[i].poly1    = scratch + i;
                args1[i].poly2    = v;
                args1[i].poly2inv = vinv;

                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                            _nmod_poly_precompute_matrix_worker, &args1[i]);
            }

            for (i = 1; i < c1; i++)
                thread_pool_wait(global_thread_pool, threads[i - 1]);

            nmod_poly_rem(tmp, H + num_threads, v);

            for (i = 0; i < c1; i++)
            {
                nmod_poly_fit_length(H + num_threads + 1 + i, v->length - 1);
                _nmod_poly_set_length(H + num_threads + 1 + i, v->length - 1);
                flint_mpn_zero(H[num_threads + 1 + i].coeffs, v->length - 1);

                args2[i].A        = HH + i;
                args2[i].res      = H + num_threads + i + 1;
                args2[i].poly1    = tmp;
                args2[i].poly3    = v;
                args2[i].poly3inv = vinv;
            }

            for (i = 1; i < c1; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                  _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker,
                                                                    &args2[i]);
            }

            _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args2[0]);

            _nmod_poly_normalise(H + num_threads + 1);

            for (i = 1; i < c1; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _nmod_poly_normalise(H + num_threads + i + 1);
            }

            for (i = 0; i < c1; i++)
            {
                nmod_poly_fit_length(I + num_threads + i + 1, v->length - 1);
                _nmod_poly_set_length(I + num_threads + i + 1, v->length - 1);
                flint_mpn_zero(I[num_threads + i + 1].coeffs, v->length - 1);

                args3[i].baby = h;
                args3[i].H    = H + num_threads + i + 1;
                args3[i].m    = l;
                args3[i].res  = I + num_threads + i + 1;
                args3[i].v    = v;
                args3[i].vinv = vinv;
                args3[i].tmp = _nmod_vec_init(v->length - 1);
            }

            for (i = 1; i < c1; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                               _nmod_poly_interval_poly_worker, &args3[i]);
            }

            _nmod_poly_interval_poly_worker(&args3[0]);

            _nmod_poly_normalise(I + num_threads + 1);

            for (i = 1; i < c1; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _nmod_poly_normalise(I + num_threads + i + 1);
            }

            for (i = 0; i < c1; i++)
            {
               _nmod_vec_clear(args3[i].tmp);
            }

            nmod_poly_one(II);

            for (i = 0; i < c1; i++)
                nmod_poly_mulmod_preinv(II, II, I + num_threads + i + 1, v,
                                                                         vinv);

            nmod_poly_gcd(II, v, II);

            if (II->length > 1)
            {
                nmod_poly_remove(v, II);
                nmod_poly_reverse(vinv, v, v->length);
                nmod_poly_inv_series_newton(vinv, vinv, v->length);

                for (i = 0; i < c1; i++)
                {
                    nmod_poly_gcd(I + num_threads + i + 1,
                                                  I + num_threads + i + 1, II);
                    if (I[num_threads + i + 1].length > 1)
                        nmod_poly_remove(II, I + num_threads + i + 1);
                }
            } else
            {
                for (i = 0; i < c1; i++)
                    nmod_poly_one(I + num_threads + i + 1);
            }

            d = d + c1*l;

            if (v->length - 1 < 2*d)
                break;
        } else if (j*(num_threads + 1) < m)
        {
            for (c2 = 0, i = 0; i < num_threads + 1 && j*(num_threads + 1) + i < m; i++, c2++)
            {
                if (HH[i].c > v->length - 1)
                {
                    _nmod_poly_reduce_matrix_mod_poly(HH + num_threads + 1,
                                                                    HH + i, v);

                    nmod_mat_clear(HH + i);
                    nmod_mat_init_set(HH + i, HH + num_threads + 1);
                    nmod_mat_clear(HH + num_threads + 1);
                }
            }

            nmod_poly_rem(tmp, H + j *(num_threads + 1) - 1, v);

            for (i = 0; i < c2; i++)
            {
                nmod_poly_fit_length(H + j*(num_threads + 1) + i,
                                                                v->length - 1);
                _nmod_poly_set_length(H + j*(num_threads + 1) + i,
                                                                v->length - 1);
                flint_mpn_zero(H[j*(num_threads + 1) + i].coeffs,
                                                                v->length - 1);

                args2[i].A        = HH + i;
                args2[i].res      = H + j*(num_threads + 1) + i;
                args2[i].poly1    = tmp;
                args2[i].poly3    = v;
                args2[i].poly3inv = vinv;
            }

            for (i = 1; i < c2; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                    _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker,
                                                                    &args2[i]);
            }

            _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args2[0]);
            _nmod_poly_normalise(H + j*(num_threads + 1));

            for (i = 1; i < c2; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _nmod_poly_normalise(H + j*(num_threads + 1) + i);
            }

            for (i = 0; i < c2; i++)
            {
                nmod_poly_fit_length(I + j*(num_threads + 1) + i,
                                                                v->length - 1);
                _nmod_poly_set_length(I + j*(num_threads + 1) + i,
                                                                v->length - 1);
                flint_mpn_zero(I[j*(num_threads + 1) + i].coeffs,
                                                                v->length - 1);
                args3[i].baby = h;
                args3[i].H    = H + j*(num_threads + 1) + i;
                args3[i].m    = l;
                args3[i].res  = I + j*(num_threads + 1) + i;
                args3[i].v    = v;
                args3[i].vinv = vinv;
                args3[i].tmp = _nmod_vec_init(v->length - 1);
            }

            for (i = 1; i < c2; i++)
            {
                thread_pool_wake(global_thread_pool, threads[i - 1], 0,
                               _nmod_poly_interval_poly_worker, &args3[i]);
            }

            _nmod_poly_interval_poly_worker(&args3[0]);
            _nmod_poly_normalise(I + j*(num_threads + 1));

            for (i = 1; i < c2; i++)
            {
                thread_pool_wait(global_thread_pool, threads[i - 1]);
                _nmod_poly_normalise(I + j*(num_threads + 1) + i);
            }

            for (i = 0; i < c2; i++)
            {
               _nmod_vec_clear(args3[i].tmp);
            }

            nmod_poly_one(II);

            for (i = 0; i < c2; i++)
                nmod_poly_mulmod_preinv(II, II, I + j*(num_threads + 1) + i,
                                                                      v, vinv);

            nmod_poly_gcd(II, v, II);

            if (II->length > 1)
            {
                nmod_poly_remove(v, II);

                nmod_poly_reverse(vinv, v, v->length);
                nmod_poly_inv_series_newton(vinv, vinv, v->length);

                for (i = 0; i < c2; i++)
                {
                    nmod_poly_gcd(I + j*(num_threads + 1) + i,
                                              I + j*(num_threads + 1) + i, II);

                    if (I[j*(num_threads + 1) + i].length > 1)
                        nmod_poly_remove(II, I + j*(num_threads + 1) + i);
                }
            } else
            {
                for (i = 0; i < c2; i++)
                    nmod_poly_one(I + j*(num_threads + 1) + i);
            }

            d = d + c2*l;

            if (v->length - 1 < 2*d)
                break;
        }
    }

    if (v->length > 1)
    {
        nmod_poly_factor_insert(res, v, 1);
        (*degs)[index++] = v->length - 1;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j].length - 1 > (j + 1)*l || j == 0)
        {
            nmod_poly_set(g, I + j);

            for (i = l - 1; i >= 0 && g->length > 1; i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                nmod_poly_sub(tmp, H + j, h + i);
                nmod_poly_gcd(f, g, tmp);

                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    nmod_poly_make_monic(f, f);

                    nmod_poly_factor_insert(res, f, 1);

                    (*degs)[index++] = l*(j + 1) - i;

                    nmod_poly_remove(g, f);
                }
            }
        } else if (I[j].length > 1)
        {
            nmod_poly_make_monic(I + j, I + j);

            nmod_poly_factor_insert(res, I + j, 1);

            (*degs)[index++] = I[j].length-1;
        }
    }

    flint_give_back_threads(threads, num_threads);

    /* cleanup */
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(v);
    nmod_poly_clear(vinv);
    nmod_poly_clear(tmp);
    nmod_poly_clear(II);

    for (i = 0; i < 2*m + l + num_threads + 2; i++)
        nmod_poly_clear(h + i);

    for (i = 0; i < c1; i++)
        nmod_mat_clear(HH + i);

    flint_free(h);
    flint_free(HH);
    flint_free(args1);
    flint_free(args2);
    flint_free(args3);
}
