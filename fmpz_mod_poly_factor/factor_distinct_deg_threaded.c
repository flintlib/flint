/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <math.h>
#include <pthread.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "fmpz_mod_poly.h"

void *
_fmpz_mod_poly_interval_poly_worker(void* arg_ptr)
{
    fmpz_mod_poly_interval_poly_arg_t arg =
                               *((fmpz_mod_poly_interval_poly_arg_t *) arg_ptr);
    slong k;
    fmpz * tmp;
    fmpz_t invV;
    fmpz_init(invV);
    tmp = _fmpz_vec_init(arg.v.length - 1);

    fmpz_invmod(invV, (arg.v.coeffs + arg.v.length - 1), &arg.v.p);

    fmpz_set_ui(arg.res.coeffs, UWORD(1));

    for (k = arg.m - 1; k >= 0; k--)
    {
        _fmpz_vec_zero(tmp, arg.v.length - 1);
        if (arg.baby[k].length < arg.v.length)
          _fmpz_vec_set(tmp, arg.baby[k].coeffs, arg.baby[k].length);
        else
          _fmpz_mod_poly_rem(tmp, arg.baby[k].coeffs, arg.baby[k].length,
                             arg.v.coeffs , arg.v.length, invV, &arg.v.p);

        _fmpz_mod_poly_sub(tmp, arg.H.coeffs, arg.H.length, tmp,
                           arg.v.length - 1, &arg.v.p);

        _fmpz_mod_poly_mulmod_preinv(arg.res.coeffs, tmp, arg.v.length - 1,
                                     arg.res.coeffs, arg.v.length - 1,
                                     arg.v.coeffs, arg.v.length,
                                     arg.vinv.coeffs, arg.vinv.length,
                                     &arg.v.p);
    }

    _fmpz_vec_clear(tmp, arg.v.length - 1);
    fmpz_clear(invV);
    flint_cleanup();
    return NULL;
}

void
fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res,
                                const fmpz_mod_poly_t poly, slong * const *degs)
{
    fmpz_mod_poly_t f, g, v, vinv, tmp, II;
    fmpz_mod_poly_t *h, *H, *I, *scratch;
    slong i, j, k, l, m, n, index, d, c1 = 1, c2;
    slong num_threads = flint_get_num_threads();
    fmpz_t p;
    fmpz_mat_t * HH;
    double beta;
    pthread_t *threads;
    fmpz_mod_poly_matrix_precompute_arg_t * args1;
    fmpz_mod_poly_compose_mod_precomp_preinv_arg_t * args2;
    fmpz_mod_poly_interval_poly_arg_t * args3;

    fmpz_init(p);
    fmpz_set(p, &poly->p);
    fmpz_mod_poly_init(v, p);

    fmpz_mod_poly_make_monic(v, poly);
    
    n = fmpz_mod_poly_degree(poly);
    if (n == 1)
    {
        fmpz_mod_poly_factor_insert(res, v, 1);
        (*degs)[0] = 1;
        fmpz_mod_poly_clear(v);
        return;
    }
    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */
    fmpz_mod_poly_init(f, p);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_init(vinv, p);
    fmpz_mod_poly_init(tmp, p);
    fmpz_mod_poly_init(II, p);

    if (!(h = flint_malloc((2 * m + l + 1+ num_threads)
                           * sizeof(fmpz_mod_poly_struct))))
    {
        flint_printf("Exception (fmpz_mod_poly_factor_distinct_deg):\n");
        flint_printf("Not enough memory.\n");
        flint_abort();
    }
    H = h + (l + 1);
    I = H + m;
    scratch = I + m;
    fmpz_mod_poly_init(h[0], p);
    fmpz_mod_poly_init(h[1], p);
    for (i = 0; i < m; i++)
    {
        fmpz_mod_poly_init(H[i], p);
        fmpz_mod_poly_init(I[i], p);
    }
    for (i = 0; i < num_threads; i++)
        fmpz_mod_poly_init(scratch[i], p);

    HH      = flint_malloc(sizeof(fmpz_mat_t) * (num_threads + 1));
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args1   = flint_malloc(num_threads *
                           sizeof(fmpz_mod_poly_matrix_precompute_arg_t));
    args2   = flint_malloc(num_threads *
                        sizeof(fmpz_mod_poly_compose_mod_precomp_preinv_arg_t));
    args3   = flint_malloc(num_threads *
                           sizeof(fmpz_mod_poly_interval_poly_arg_t));

    fmpz_mod_poly_reverse(vinv, v, v->length);
    fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length);

    /* compute baby steps: h[i]=x^{p^i}mod v */
    fmpz_mod_poly_set_coeff_ui(h[0], 1, 1);
    fmpz_mod_poly_powmod_x_fmpz_preinv(h[1], p, v, vinv);
    if (fmpz_sizeinbase(p, 2) > ((n_sqrt(v->length - 1) + 1) * 3) / 4)
    {
        for (i = 1; i < FLINT_BIT_COUNT(l); i++)
            fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(*(h + 1 +
                                                             (1 << (i - 1))),
                                                             *(h + 1),
                                                             (1 << (i - 1)),
                                                             (1 << (i - 1)), v,
                                                             vinv);
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(*(h + 1 +
                                                         (1 << (i - 1))),
                                                         *(h + 1),
                                                         (1 << (i - 1)),
                                                         l - (1 << (i - 1)), v,
                                                         vinv);
    }
    else
    {
        for (i = 2; i < l + 1; i++)
        {
            fmpz_mod_poly_init(h[i], p);
            fmpz_mod_poly_powmod_fmpz_binexp_preinv(h[i], h[i - 1], p,
                                                    v, vinv);
        }
    }

    /* compute coarse distinct-degree factorisation */
    index = 0;
    fmpz_mod_poly_set(H[0], h[l]);
    fmpz_mat_init(HH[0], n_sqrt(v->length - 1) + 1, v->length - 1);
    fmpz_mod_poly_precompute_matrix(HH[0], H[0], v, vinv);

    d = 1;
    for (j = 0; j < m / num_threads + 1; j++)
    {
        if (j == 0)
        {
            for (i = 0; i < num_threads; i++)
            {
                if (i > 0 && I[i - 1]->length > 1)
                {
                    _fmpz_mod_poly_reduce_matrix_mod_poly(HH[num_threads],
                                                          HH[0], v);
                    fmpz_mat_clear(HH[0]);
                    fmpz_mat_init_set(HH[0], HH[num_threads]);
                    fmpz_mat_clear(HH[num_threads]);
                    fmpz_mod_poly_rem(tmp, H[i - 1], v);
                    fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[i],
                                                                     tmp, HH[0],
                                                                     v, vinv);
                }
                else if (i > 0)
                    fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[i],
                                                             H[i - 1], HH[0], v,
                                                             vinv);

                /* compute interval polynomials */
                fmpz_mod_poly_set_coeff_ui(I[i], 0, 1);
                for (k = l - 1; (k >= 0) && (2 * d <= v->length - 1); k--, d++)
                {
                    fmpz_mod_poly_rem(tmp, h[k], v);
                    fmpz_mod_poly_sub(tmp, H[i], tmp);
                    fmpz_mod_poly_mulmod_preinv(I[i], tmp, I[i], v, vinv);
                }

                /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
                /* F_j is stored on the place of I_j */
                fmpz_mod_poly_gcd(I[i], v, I[i]);
                if (I[i]->length > 1)
                {
                    fmpz_mod_poly_remove(v, I[i]);
                    fmpz_mod_poly_reverse(vinv, v, v->length);
                    fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length);
                }
                if (v->length - 1 < 2 * d)
                    break;
            }
            if (v->length - 1 < 2 * d)
                break;
        }
        else if (j == 1 && num_threads < m)
        {
            if (I[num_threads - 1]->length > 1)
            {
                _fmpz_mod_poly_reduce_matrix_mod_poly(HH[num_threads], HH[0],
                                                      v);
                fmpz_mat_clear(HH[0]);
                fmpz_mat_init_set(HH[0], HH[num_threads]);
                fmpz_mat_clear(HH[num_threads]);
            }

            c1 = 1;
            for (i = 1; i < num_threads && i + num_threads < m; i++, c1++)
            {
                fmpz_mat_init(HH[i], n_sqrt(v->length - 1) + 1, v->length - 1);
                fmpz_mod_poly_rem(scratch[i], H[i], v);
                if (scratch[i]->length < v->length - 1)
                {
                    fmpz_mod_poly_fit_length(scratch[i], v->length - 1);
                    _fmpz_vec_zero(scratch[i]->coeffs + scratch[i]->length,
                                   v->length - 1 - scratch[i]->length);
                    _fmpz_mod_poly_set_length(scratch[i], v->length - 1);
                }
                args1[i].A        = *HH[i];
                args1[i].poly1    = *scratch[i];
                args1[i].poly2    = *v;
                args1[i].poly2inv = *vinv;

                pthread_create(&threads[i], NULL,
                            _fmpz_mod_poly_precompute_matrix_worker, &args1[i]);
            }
            for (i = 1; i < c1; i++)
                pthread_join(threads[i], NULL);

            fmpz_mod_poly_rem(tmp, H[num_threads - 1], v);
            for (i = 0; i < c1; i++)
            {
                fmpz_mod_poly_fit_length(H[num_threads + i], v->length - 1);
                _fmpz_mod_poly_set_length(H[num_threads + i], v->length - 1);
                args2[i].A        = *HH[i];
                args2[i].res      = *H[num_threads + i];
                args2[i].poly1    = *tmp;
                args2[i].poly3    = *v;
                args2[i].poly3inv = *vinv;

                pthread_create(&threads[i], NULL,
        _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args2[i]);
            }
            for (i = 0; i < c1; i++)
            {
                pthread_join(threads[i], NULL);
                _fmpz_mod_poly_normalise(H[num_threads + i]);
            }

            for (i = 0; i < c1; i++)
            {
                fmpz_mod_poly_fit_length(I[num_threads + i], v->length - 1);
                _fmpz_mod_poly_set_length(I[num_threads + i], v->length - 1);
                _fmpz_vec_zero(I[num_threads + i]->coeffs, v->length - 1);
                args3[i].baby = *h;
                args3[i].H    = *H[num_threads + i];
                args3[i].m    = l;
                args3[i].res  = *I[num_threads + i];
                args3[i].v    = *v;
                args3[i].vinv = *vinv;

                pthread_create(&threads[i], NULL,
                               _fmpz_mod_poly_interval_poly_worker, &args3[i]);
            }

            for (i = 0; i < c1; i++)
            {
                pthread_join(threads[i], NULL);
                _fmpz_mod_poly_normalise(I[num_threads + i]);
            }

            fmpz_mod_poly_set_ui(II, UWORD(1));

            for (i = 0; i < c1; i++)
                fmpz_mod_poly_mulmod_preinv(II, II, I[num_threads + i], v,
                                            vinv);

            fmpz_mod_poly_gcd(II, v, II);
            if (II->length > 1)
            {
                fmpz_mod_poly_remove(v, II);
                fmpz_mod_poly_reverse(vinv, v, v->length);
                fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length);
                for (i = 0; i < c1; i++)
                {
                    fmpz_mod_poly_gcd(I[num_threads + i], I[num_threads + i],
                                      II);
                    if (I[num_threads + i]->length > 1)
                        fmpz_mod_poly_remove(II, I[num_threads + i]);
                }
            }
            else
            {
                for (i = 0; i < c1; i++)
                    fmpz_mod_poly_set_ui(I[num_threads + i], UWORD(1));
            }
            d = d + c1 * l;
            if (v->length-1 < 2 * d)
                break;
        }
        else if (j*num_threads < m)
        {
            c2 = 0;
            for (i = 0; i < num_threads && j*num_threads + i < m; i++, c2++)
            {
                if (HH[i] -> c > v -> length - 1)
                {
                    _fmpz_mod_poly_reduce_matrix_mod_poly(HH[num_threads],
                                                          HH[i], v);
                    fmpz_mat_clear(HH[i]);
                    fmpz_mat_init_set(HH[i], HH[num_threads]);
                    fmpz_mat_clear(HH[num_threads]);
                }
            }

            fmpz_mod_poly_rem(tmp, H[j * num_threads - 1], v);
            for (i = 0; i < c2; i++)
            {
                fmpz_mod_poly_fit_length(H[j * num_threads + i], v->length - 1);
                _fmpz_mod_poly_set_length(H[j * num_threads + i],
                                          v->length - 1);
                args2[i].A        = *HH[i];
                args2[i].res      = *H[j * num_threads + i];
                args2[i].poly1    = *tmp;
                args2[i].poly3    = *v;
                args2[i].poly3inv = *vinv;

                pthread_create(&threads[i], NULL,
        _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args2[i]);
            }
            for (i = 0; i < c2; i++)
            {
                pthread_join(threads[i], NULL);
                _fmpz_mod_poly_normalise(H[j * num_threads + i]);
            }

            for (i = 0; i < c2; i++)
            {
                fmpz_mod_poly_fit_length(I[j * num_threads + i], v->length - 1);
                _fmpz_mod_poly_set_length(I[j * num_threads + i],
                                          v->length - 1);
                _fmpz_vec_zero(I[j * num_threads + i]->coeffs, v->length - 1);
                args3[i].baby = *h;
                args3[i].H    = *H[j * num_threads + i];
                args3[i].m    = l;
                args3[i].res  = *I[j * num_threads + i];
                args3[i].v    = *v;
                args3[i].vinv = *vinv;

                pthread_create(&threads[i], NULL,
                               _fmpz_mod_poly_interval_poly_worker, &args3[i]);
            }

            for (i = 0; i < c2; i++)
            {
                pthread_join(threads[i], NULL);
                _fmpz_mod_poly_normalise(I[j * num_threads + i]);
            }

            fmpz_mod_poly_set_ui(II, UWORD(1));

            for (i = 0; i < c2; i++)
                fmpz_mod_poly_mulmod_preinv(II, II, I[j * num_threads + i], v,
                                            vinv);

            fmpz_mod_poly_gcd(II, v, II);
            if (II->length > 1)
            {
                fmpz_mod_poly_remove(v, II);
                fmpz_mod_poly_reverse(vinv, v, v->length);
                fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length);
                for (i = 0; i < c2; i++)
                {
                    fmpz_mod_poly_gcd(I[j * num_threads + i],
                                      I[j * num_threads + i], II);
                    if (I[j * num_threads + i]->length > 1)
                        fmpz_mod_poly_remove(II, I[j * num_threads + i]);
                }
            }
            else
            {
                for (i = 0; i < c2; i++)
                    fmpz_mod_poly_set_ui(I[j * num_threads + i], UWORD(1));
            }
            d = d + c2 * l;
            if (v->length - 1 < 2 * d)
                break;
        }
    }
    if (v->length > 1)
    {
        fmpz_mod_poly_factor_insert(res, v, 1);
        (*degs)[index++] = v->length - 1;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length - 1 > (j + 1)*l || j == 0)
        {
            fmpz_mod_poly_set(g, I[j]);
            for (i = l - 1; i >= 0 && (g->length > 1); i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                fmpz_mod_poly_sub(tmp, H[j], h[i]);
                fmpz_mod_poly_gcd(f, g, tmp);
                if (f->length > 1)
                {
                    fmpz_mod_poly_make_monic(f, f);
                    fmpz_mod_poly_factor_insert(res, f, 1);
                    (*degs)[index++] = l * (j + 1) - i;

                    fmpz_mod_poly_remove(g, f);
                }
            }
        }
        else if (I[j]->length > 1)
        {
            fmpz_mod_poly_make_monic(I[j], I[j]);
            fmpz_mod_poly_factor_insert(res, I[j], 1);
            (*degs)[index++] = I[j]->length - 1;
        }
    }

    /* cleanup */
    fmpz_clear(p);
    fmpz_mod_poly_clear(f);
    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_clear(II);
    fmpz_mod_poly_clear(vinv);
    fmpz_mod_poly_clear(tmp);

    for (i = 0; i < l + 1; i++)
        fmpz_mod_poly_clear(h[i]);
    for (i = 0; i < m; i++)
    {
        fmpz_mod_poly_clear(H[i]);
        fmpz_mod_poly_clear(I[i]);
    }
    for (i = 0; i < num_threads; i++)
        fmpz_mod_poly_clear(scratch[i]);
    for (i = 0; i < c1; i++)
        fmpz_mat_clear(HH[i]);

    flint_free(h);
    flint_free(HH);
    flint_free(args1);
    flint_free(args2);
    flint_free(args3);
    flint_free(threads);
}
