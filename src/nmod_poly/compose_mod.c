/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

void
_nmod_poly_compose_mod(mp_ptr res,
    mp_srcptr f, slong lenf, mp_srcptr g, mp_srcptr h, slong lenh, nmod_t mod)
{
    if (lenh < 8 || lenf >= lenh)
        _nmod_poly_compose_mod_horner(res, f, lenf, g, h, lenh, mod);
    else
        _nmod_poly_compose_mod_brent_kung(res, f, lenf, g, h, lenh, mod);
}

void
nmod_poly_compose_mod(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod). Division by zero.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}

void
_nmod_poly_compose_mod_brent_kung(mp_ptr res, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2,
                            mp_srcptr poly3, slong len3, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }

    if (len3 == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
        return;
    }

    m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    A->rows[0][0] = UWORD(1);
    _nmod_vec_set(A->rows[1], poly2, n);
    for (i = 2; i < m; i++)
        _nmod_poly_mulmod(A->rows[i], A->rows[i-1],
            n, poly2, n, poly3, len3, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod(h, A->rows[m - 1], n, poly2, n, poly3, len3, mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod(t, res, n, h, n, poly3, len3, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_brent_kung). The degree of the \n"
               "first polynomial must be smaller than that of the modulus.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod_brent_kung(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_brent_kung(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}

void
_nmod_poly_reduce_matrix_mod_poly(nmod_mat_t A, const nmod_mat_t B,
                                   const nmod_poly_t f)
{
    mp_ptr tmp1;
    slong n = f->length - 1;
    slong i, m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, f->mod.n);

    tmp1 = _nmod_vec_init(B->c - f->length + 1);

    A->rows[0][0] = 1;

    for (i = 1; i < m; i++)
        _nmod_poly_divrem(tmp1, A->rows[i], B->rows[i], B->c, f->coeffs,
                                                            f->length, f->mod);

    _nmod_vec_clear(tmp1);
}

void
_nmod_poly_precompute_matrix(nmod_mat_t A, mp_srcptr poly1, mp_srcptr poly2,
                     slong len2, mp_srcptr poly2inv, slong len2inv, nmod_t mod)
{
    /* Set rows of A to powers of poly1 */
    slong n, m;

    n = len2 - 1;

    m = n_sqrt(n) + 1;

    _nmod_poly_powers_mod_preinv_naive(A->rows, poly1, n, m,
                                          poly2, len2, poly2inv, len2inv, mod);
}

void
nmod_poly_precompute_matrix(nmod_mat_t A, const nmod_poly_t poly1,
                            const nmod_poly_t poly2, const nmod_poly_t poly2inv)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = len2 - 1;
    slong m = n_sqrt(len) + 1;

    mp_ptr ptr1;

    if (len2 == 0)
    {
        flint_printf("Exception (nmod_poly_precompute_matrix). Division by zero.\n");
        flint_abort();
    }

    if (A->r != m || A->c != len)
    {
        flint_printf("Exception (nmod_poly_precompute_matrix). Wrong dimensions.\n");
        flint_abort();
    }

    if (len2 == 1)
    {
        nmod_mat_zero(A);
        return;
    }

    ptr1 = _nmod_vec_init(len);

    if (len1 <= len)
    {
        flint_mpn_copyi(ptr1, poly1->coeffs, len1);
        flint_mpn_zero(ptr1 + len1, len - len1);
    } else
    {
        _nmod_poly_rem(ptr1, poly1->coeffs, len1, poly2->coeffs, len2, A->mod);
    }

    _nmod_poly_precompute_matrix(A, ptr1, poly2->coeffs,
                             len2, poly2inv->coeffs, poly2inv->length, A->mod);

    _nmod_vec_clear(ptr1);
}

void
_nmod_poly_compose_mod_brent_kung_precomp_preinv(mp_ptr res, mp_srcptr poly1,
                  slong len1, const nmod_mat_t A, mp_srcptr poly3, slong len3,
                                 mp_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nmod_mat_t B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }

    if (len3 == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, A->rows[1][0], mod);
        return;
    }

    m = n_sqrt(n) + 1;

    /* TODO check A*/

    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1/m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1%m);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n,
                                           poly3, len3, poly3inv, len3inv,mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                                       poly3inv, len3inv, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). The degree of the \n"
               "first polynomial must be smaller than that of the modulus.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);

        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);

        return;
    }

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        nmod_poly_t tmp;

        nmod_poly_init_mod(tmp, res->mod);

        nmod_poly_compose_mod_brent_kung_precomp_preinv(tmp, poly1, A,
                                                              poly3, poly3inv);

        nmod_poly_swap(tmp, res);

        nmod_poly_clear(tmp);

        return;
    }

    nmod_poly_fit_length(res, len);

    _nmod_poly_compose_mod_brent_kung_precomp_preinv(res->coeffs,
                            poly1->coeffs, len1, A, poly3->coeffs, len3,
                                 poly3inv->coeffs, poly3inv->length, res->mod);

    res->length = len;

    _nmod_poly_normalise(res);
}

void
_nmod_poly_compose_mod_brent_kung_preinv(mp_ptr res, mp_srcptr poly1,
                  slong len1, mp_srcptr poly2, mp_srcptr poly3, slong len3,
                                 mp_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }

    if (len3 == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
        return;
    }

    m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1/m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1%m);

    /* Set rows of A to powers of poly2 */
    _nmod_poly_powers_mod_preinv_naive(A->rows, poly2, n,
                                       m, poly3, len3, poly3inv, len3inv, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, poly2, n,
                                           poly3, len3, poly3inv, len3inv,mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_preinv(t, res, n, h, n, poly3, len3,
                                                       poly3inv, len3inv, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                           const nmod_poly_t poly3, const nmod_poly_t poly3inv)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_preinv). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_preinv). The degree of the \n"
               "first polynomial must be smaller than that of the modulus.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1 || res == poly3inv)
    {
        nmod_poly_t tmp;
        nmod_poly_init_mod(tmp, res->mod);
        nmod_poly_compose_mod_brent_kung_preinv(tmp, poly1, poly2,
                                                              poly3, poly3inv);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    } else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);

    _nmod_poly_compose_mod_brent_kung_preinv(res->coeffs,
                      poly1->coeffs, len1, ptr2, poly3->coeffs, len3,
                                 poly3inv->coeffs, poly3inv->length, res->mod);

    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}

void
_nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                const nmod_poly_struct * polys, slong lenpolys, slong l,
                         mp_srcptr g, slong glen, mp_srcptr poly, slong len,
                                   mp_srcptr polyinv, slong leninv, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, j, k, n, m, len2 = l, len1;

    n = len - 1;

    m = n_sqrt(n*len2) + 1;

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    k = len/m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k*len2, m, mod.n);
    nmod_mat_init(C, k*len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = (polys + j)->length;

        for (i = 0; i < len1/m; i++)
            _nmod_vec_set(B->rows[i + j*k], (polys + j)->coeffs + i*m, m);

        _nmod_vec_set(B->rows[i + j*k], (polys + j)->coeffs + i*m, len1%m);
    }

    /* Set rows of A to powers of last element of polys */
    _nmod_poly_powers_mod_preinv_naive(A->rows, g, glen,
                                           m, poly, len, polyinv, leninv, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        h[0] = n_mulmod2_preinv(A->rows[m - 1][0],
                                               A->rows[1][0], mod.n, mod.ninv);
    } else
    {
        _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                                                    len, polyinv, leninv, mod);
    }

    for (j = 0; j < len2; j++)
    {
        _nmod_vec_set((res + j)->coeffs, C->rows[(j + 1)*k - 1], n);

        if (n == 1)
        {
            for (i = 2; i <= k; i++)
            {
                t[0] = n_mulmod2_preinv(res[j].coeffs[0],
                                                        h[0], mod.n, mod.ninv);
                res[j].coeffs[0] = n_addmod(t[0],
                                             C->rows[(j + 1)*k - i][0], mod.n);
            }
        } else
        {
            for (i = 2; i <= k; i++)
            {
                _nmod_poly_mulmod_preinv(t, res[j].coeffs,
                                     n, h, n, poly, len, polyinv, leninv, mod);
                _nmod_poly_add(res[j].coeffs, t, n,
                                               C->rows[(j + 1)*k - i], n, mod);
            }
        }
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                     const nmod_poly_struct * polys, slong len1, slong n,
        const nmod_poly_t g, const nmod_poly_t poly, const nmod_poly_t polyinv)
{
    slong len2 = poly->length;
    slong len3, i;

    for (i = 0; i < len1; i++)
    {
        len3 = (polys + i)->length;

        if (len3 >= len2)
        {
            flint_printf("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv)."
                 "The degree of the first polynomial must be smaller than that of the "
                 " modulus\n");
            flint_abort();
        }
    }

    if (n > len1)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv)."
                                     "n is larger than the length of polys\n");
        flint_abort();
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n,
                g->coeffs, g->length, poly->coeffs, len2, polyinv->coeffs,
                                                   polyinv->length, poly->mod);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}

typedef struct
{
    nmod_poly_struct * res;
    nmod_mat_struct * C;
    mp_srcptr h;
    mp_srcptr poly;
    mp_srcptr polyinv;
    nmod_t p;
    mp_ptr t;
    volatile slong * j;
    slong k;
    slong m;
    slong len;
    slong leninv;
    slong len2;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
} compose_vec_arg_t;

void
_nmod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr)
{
    compose_vec_arg_t arg = *((compose_vec_arg_t *) arg_ptr);
    slong i, j, k = arg.k, n = arg.len - 1;
    slong len = arg.len, leninv = arg.leninv;
    mp_ptr t = arg.t;
    mp_srcptr h = arg.h;
    mp_srcptr poly = arg.poly;
    mp_srcptr polyinv = arg.polyinv;
    nmod_poly_struct * res = arg.res;
    nmod_mat_struct * C = arg.C;
    nmod_t p = arg.p;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        j = *arg.j;
        *arg.j = j + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (j >= arg.len2)
            return;

        _nmod_vec_set(res[j].coeffs, C->rows[(j + 1)*k - 1], n);

        if (n == 1) /* special case, constant polynomials */
        {
            for (i = 2; i <= k; i++)
            {
                t[0] = n_mulmod2_preinv(res[j].coeffs[0], h[0], p.n, p.ninv);
                res[j].coeffs[0] = n_addmod(t[0],
                                               C->rows[(j + 1)*k - i][0], p.n);
            }
        } else
        {
            for (i = 2; i <= k; i++)
            {
                _nmod_poly_mulmod_preinv(t, res[j].coeffs, n, h, n, poly,
                                                      len, polyinv, leninv, p);
                _nmod_poly_add(res[j].coeffs, t, n,
                                                 C->rows[(j + 1)*k - i], n, p);
            }
        }
    }
}

void
_nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
                                             const nmod_poly_struct * polys,
                                             slong lenpolys, slong l,
                                             mp_srcptr g, slong glen,
                                             mp_srcptr poly, slong len,
                                             mp_srcptr polyinv, slong leninv,
                                             nmod_t mod,
                                             thread_pool_handle * threads,
                                             slong num_threads)
{
    nmod_mat_t A, B, C;
    slong i, j, n, m, k, len2 = l, len1, shared_j = 0;
    mp_ptr h;
    compose_vec_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    n = len - 1;

    m = n_sqrt(n*len2) + 1;

    h = _nmod_vec_init(n);

    k = len/m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k*len2, m, mod.n);
    nmod_mat_init(C, k*len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = polys[j].length;

        for (i = 0; i < len1 / m; i++)
            _nmod_vec_set(B->rows[i + j * k], polys[j].coeffs + i * m, m);

        _nmod_vec_set(B->rows[i + j * k], polys[j].coeffs + i * m,
                      len1 % m);
    }

    /* Set rows of A to powers of g */
    _nmod_poly_powers_mod_preinv_threaded_pool(A->rows, g, glen,
                         m, poly, len, polyinv, leninv, mod, threads, num_threads);

    _nmod_mat_mul_classical_threaded_pool_op(C, NULL, B, A, 0,
                                                         threads, num_threads);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        h[0] = n_mulmod2_preinv(A->rows[m - 1][0],
                                               A->rows[1][0], mod.n, mod.ninv);
    } else
    {

        _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                             len, polyinv, leninv, mod);
    }

    args = (compose_vec_arg_t *)
                   flint_malloc(sizeof(compose_vec_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].res     = res;
        args[i].C       = C;
        args[i].h       = h;
        args[i].k       = k;
        args[i].m       = m;
        args[i].j       = &shared_j;
        args[i].poly    = poly;
        args[i].t       = _nmod_vec_init(len);
        args[i].len     = len;
        args[i].polyinv = polyinv;
        args[i].leninv  = leninv;
        args[i].p       = mod;
        args[i].len2    = len2;
#if FLINT_USES_PTHREAD
        args[i].mutex   = &mutex;
#endif
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _nmod_poly_compose_mod_brent_kung_vec_preinv_worker, &args[i]);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    for (i = 0; i < num_threads + 1; i++)
       _nmod_vec_clear(args[i].t);

    flint_free(args);

    _nmod_vec_clear(h);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t g,
                                            const nmod_poly_t poly,
                                            const nmod_poly_t polyinv,
                                            thread_pool_handle * threads,
                                            slong num_threads)
{
    slong len2 = poly->length;
    slong i;

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys, len1, n,
                                                          g->coeffs, g->length,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          poly->mod,
                                                          threads,
                                                          num_threads);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t g,
                                            const nmod_poly_t poly,
                                            const nmod_poly_t polyinv)
{
    slong i, len2 = poly->length;
    thread_pool_handle * threads;
    slong num_threads;

    for (i = 0; i < len1; i++)
    {
        if (polys[i].length >= len2)
        {
            flint_printf
                ("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_threaded)."
                 "The degree of the first polynomial must be smaller than that of the "
                 " modulus\n");
            flint_abort();
        }
    }

    if (n > len1)
    {
        flint_printf
            ("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_threaded)."
             "n is larger than the length of polys\n");
        flint_abort();
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys, len1, n,
                                                          g->coeffs, g->length,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          poly->mod,
                                                          threads,
                                                          num_threads);

    flint_give_back_threads(threads, num_threads);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}

void
_nmod_poly_compose_mod_horner(mp_ptr res,
    mp_srcptr f, slong lenf, mp_srcptr g, mp_srcptr h, slong lenh, nmod_t mod)
{
    slong i, len;
    mp_ptr t;

    if (lenh == 1)
        return;

    if (lenf == 1)
    {
        res[0] = f[0];
        return;
    }

    if (lenh == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(f, lenf, g[0], mod);
        return;
    }

    len = lenh - 1;
    i = lenf - 1;
    t = _nmod_vec_init(len);

    _nmod_vec_scalar_mul_nmod(res, g, len, f[i], mod);
    i--;
    if (i >= 0)
        res[0] = nmod_add(res[0], f[i], mod);

    while (i > 0)
    {
        i--;
        _nmod_poly_mulmod(t, res, len, g, len, h, lenh, mod);
        _nmod_poly_add(res, t, len, f + i, 1, mod);
    }

    _nmod_vec_clear(t);
}

void
nmod_poly_compose_mod_horner(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_horner). Division by zero.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod_horner(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_horner(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
