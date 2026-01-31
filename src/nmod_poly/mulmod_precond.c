/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly/impl.h"

/* This should match _nmod_poly_divrem_try_sparse. Todo:
   implement a sparse method in this file instead of relying on
   _nmod_poly_divrem_newton_n_preinv to have a sparse hack. */
#define MAX_NZ 6

static int
_nmod_poly_mulmod_precond_select_method(nn_srcptr d, slong dlen, slong num, nmod_t mod)
{
    int sparse;
    slong i, nz;
    slong n = dlen - 1;

    if (n == 1) return NMOD_POLY_MULMOD_PRECOND_NONE;
    if (n <= 16 && num >= 2) return NMOD_POLY_MULMOD_PRECOND_MATRIX;

    sparse = 1;
    nz = 1;
    for (i = 0; i < dlen - 1 && sparse; i++)
    {
        nz += (d[i] != 0);
        sparse = nz < MAX_NZ;
    }

    if (!sparse)
    {
        if (n <= 64 && num >= 4) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
        if (n <= 128 && num >= 8) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
        if (n <= 192 && num >= 32) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
        if (n <= 320 && num >= 64 && mod.n <= 255) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
        if (num >= 8) return NMOD_POLY_MULMOD_PRECOND_SHOUP;
    }
    else
    {
        if (n <= 32 && num >= 8) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
        if (n <= 64 && num >= 16) return NMOD_POLY_MULMOD_PRECOND_MATRIX;
    }

    return NMOD_POLY_MULMOD_PRECOND_NONE;
}

static void
_nmod_poly_mulmod_precond_init_none(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, nmod_t FLINT_UNUSED(mod))
{
    precond->method = NMOD_POLY_MULMOD_PRECOND_NONE;
    precond->a = a;
    precond->alen = alen;
    precond->n = dlen - 1;
    precond->d = d;
    precond->dinv = dinv;
    precond->lendinv = lendinv;
}

static void
_nmod_poly_mulmod_precond_init_shoup(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, nmod_t mod)
{
    slong n;
    nn_ptr arev;
    TMP_INIT;

    precond->method = NMOD_POLY_MULMOD_PRECOND_SHOUP;

    n = dlen - 1;

    precond->n = n;
    precond->a = a;
    precond->alen = alen;
    precond->d = d;

    FLINT_ASSERT(alen >= 0);
    FLINT_ASSERT(dlen >= 1);
    FLINT_ASSERT(alen < dlen);

    precond->adivd = _nmod_vec_init(n);

    TMP_START;
    arev = TMP_ALLOC(n * sizeof(ulong));

    /* Todo: could avoid zero-padding when alen < n. Also, we actually just
       need the quotient to precision n - 1, so we could compute one less
       coefficient here. */
    _nmod_poly_reverse(arev, a, alen, n);
    _nmod_poly_mullow(precond->adivd, arev, n, dinv, FLINT_MIN(n, lendinv), n, mod);
    _nmod_poly_reverse(precond->adivd, precond->adivd, n, n);

    TMP_END;
}

static void
_nmod_poly_mulmod_precond_shoup(nn_ptr res, nn_srcptr a, slong alen, nn_srcptr adivd, slong n, nn_srcptr b, slong blen, nn_srcptr d, nmod_t mod)
{
    nn_ptr T, U, V;
    TMP_INIT;

    /* needs a special case */
    if (alen + blen - 1 <= n)
    {
        if (alen >= blen)
            _nmod_poly_mul(res, a, alen, b, blen, mod);
        else
            _nmod_poly_mul(res, b, blen, a, alen, mod);
        _nmod_vec_zero(res + alen + blen - 1, n - (alen + blen - 1));
        return;
    }

    FLINT_ASSERT(alen <= n);
    FLINT_ASSERT(blen <= n);

    TMP_START;
    T = TMP_ALLOC(sizeof(ulong) * (3 * n));
    U = T + n;
    V = U + n;

    /* n - 1 high coefficients of ab/d */
    _nmod_poly_reverse(T, adivd + 1, n - 1, n - 1);
    _nmod_poly_reverse(U, b + 1, blen - 1, n - 1);
    _nmod_poly_mullow(V, T, n - 1, U, n - 1, n - 1, mod);
    _nmod_poly_reverse(V, V, n - 1, n - 1);

    /* low multiplication by d */
    _nmod_poly_mullow(U, d, n, V, n - 1, n, mod);

    if (alen >= blen)
        _nmod_poly_mullow(T, a, alen, b, blen, n, mod);
    else
        _nmod_poly_mullow(T, b, blen, a, alen, n, mod);

    _nmod_vec_sub(res, T, U, n, mod);

    TMP_END;
}

static void
_nmod_poly_mulmod_precond_precompute_matrix(nn_ptr res, nn_srcptr a, slong alen, slong n, nn_srcptr d, ulong dinv, int packing, nmod_t mod)
{
    slong i, j, k;
    ulong q[1];
    nn_srcptr p;
    nn_ptr t;
    TMP_INIT;

    FLINT_ASSERT(alen <= n);
    FLINT_ASSERT(1 <= packing && packing <= 4);

    TMP_START;
    t = TMP_ALLOC(sizeof(ulong) * (2 * n));

#define RES(ii,jj) ((res)[(ii)*n + (jj)])

    if (packing == 1)
    {
        for (j = 0; j < alen; j++)
            RES(j, 0) = t[n + j] = a[j];
        for ( ; j < n; j++)
            RES(j, 0) = t[n + j] = 0;

        for (i = 1; i < n; i++)
        {
            t[n - i] = 0;
            _nmod_poly_divrem_q0_preinv1(q, t + n - i, t + n - i, d, n + 1, dinv, mod);
            for (j = 0; j < n; j++)
                RES(j, i) = t[n - i + j];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            if (i == 0)
            {
                for (j = 0; j < alen; j++)
                    t[n + j] = a[j];
                for ( ; j < n; j++)
                    t[n + j] = 0;
            }
            else
            {
                t[n - i] = 0;
                _nmod_poly_divrem_q0_preinv1(q, t + n - i, t + n - i, d, n + 1, dinv, mod);
            }

            if (packing == 2)
            {
                for (j = k = 0; j + 1 < n; j += 2, k++)
                {
                    p = t + n - i + j;
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 2));
                }

                p = t + n - i + j;
                if (j == n - 1)
                    RES(k, i) = p[0];
            }
            else if (packing == 3)
            {
                for (j = k = 0; j + 2 < n; j += 3, k++)
                {
                    p = t + n - i + j;
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 3))
                                     | (p[2] << (2 * (FLINT_BITS / 3)));
                }

                p = t + n - i + j;
                if (j == n - 1)
                    RES(k, i) = p[0];
                else if (j == n - 2)
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 3));
            }
            else if (packing == 4)
            {
                for (j = k = 0; j + 3 < n; j += 4, k++)
                {
                    p = t + n - i + j;
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 4))
                                     | (p[2] << (2 * (FLINT_BITS / 4)))
                                     | (p[3] << (3 * (FLINT_BITS / 4)));
                }

                p = t + n - i + j;
                if (j == n - 1)
                    RES(k, i) = p[0];
                else if (j == n - 2)
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 4));
                else if (j == n - 3)
                    RES(k, i) = p[0] | (p[1] << (FLINT_BITS / 4))
                                     | (p[2] << (2 * (FLINT_BITS / 4)));
            }
        }
    }
#undef RES

    TMP_END;
}

static void
_nmod_poly_mulmod_precond_init_matrix(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, ulong dinv, nmod_t mod)
{
    slong n;
    ulong m;
    int packing;
    slong alloc;

    precond->method = NMOD_POLY_MULMOD_PRECOND_MATRIX;

    FLINT_ASSERT(alen >= 0);
    FLINT_ASSERT(dlen >= 1);
    FLINT_ASSERT(alen < dlen);

    n = dlen - 1;
    m = mod.n;

    precond->n = n;
    precond->a = a;
    precond->alen = alen;

    if (n == 0)
    {
        precond->packing = 1;
        precond->dot_params = _nmod_vec_dot_params(0, mod);
        precond->matrix = NULL;
        return;
    }

    packing = 1;
    alloc = n * n;

    if (m < (UWORD(1) << (FLINT_BITS / 2)))
    {
        ulong hi, bound;
        m = (m - 1) * (m - 1);
        umul_ppmm(hi, bound, m, n);

        if (hi == 0)
        {
            if (bound < (UWORD(1) << (FLINT_BITS / 4)))
            {
                packing = 4;
                alloc = ((n + 3) / 4) * n;
            }
            else if (bound < (UWORD(1) << (FLINT_BITS / 3)))
            {
                packing = 3;
                alloc = ((n + 2) / 3) * n;
            }
            else if (bound < (UWORD(1) << (FLINT_BITS / 2)))
            {
                packing = 2;
                alloc = ((n + 1) / 2) * n;
            }
        }
    }

    precond->packing = packing;

    if (packing == 1)
        precond->dot_params = _nmod_vec_dot_params(n, mod);

    precond->matrix = flint_malloc(sizeof(ulong) * alloc);
    _nmod_poly_mulmod_precond_precompute_matrix(precond->matrix, a, alen, n, d, dinv, packing, mod);
}

static void
_nmod_poly_mulmod_precond_matrix(nn_ptr res, nn_srcptr apre, slong len, nn_srcptr b, slong blen, int packing, dot_params_t params, nmod_t mod)
{
    slong i, j;
    ulong c, c0, c1, c2, c3;
    slong n = len;

    if (packing == 1)
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod_vec_dot(apre + i * len, b, blen, mod, params);
    }
    else if (packing == 2)
    {
        for (i = 0; i + 1 < n; i += 2)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 2) * n + j] * b[j];

            c0 = c % (UWORD(1) << (FLINT_BITS / 2));
            c1 = (c >> (FLINT_BITS / 2));

            res[i + 0] = nmod_set_ui(c0, mod);
            res[i + 1] = nmod_set_ui(c1, mod);
        }

        if (n % 2)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 2) * n + j] * b[j];

            c0 = c;
            res[i + 0] = nmod_set_ui(c0, mod);
        }
    }
    else if (packing == 3)
    {
        for (i = 0; i + 2 < n; i += 3)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 3) * n + j] * b[j];

            c0 = c % (UWORD(1) << (FLINT_BITS / 3));
            c1 = (c >> (FLINT_BITS / 3)) % (UWORD(1) << (FLINT_BITS / 3));
            c2 = (c >> (2 * (FLINT_BITS / 3)));

            res[i + 0] = nmod_set_ui(c0, mod);
            res[i + 1] = nmod_set_ui(c1, mod);
            res[i + 2] = nmod_set_ui(c2, mod);
        }

        if (n % 3)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 3) * n + j] * b[j];

            c0 = c % (UWORD(1) << (FLINT_BITS / 3));
            res[i + 0] = nmod_set_ui(c0, mod);

            if (n % 3 == 2)
            {
                c1 = c >> (FLINT_BITS / 3);
                res[i + 1] = nmod_set_ui(c1, mod);
            }
        }
    }
    else
    {
        for (i = 0; i + 3 < n; i += 4)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 4) * n + j] * b[j];

            c0 = c % (UWORD(1) << (FLINT_BITS / 4));
            c1 = (c >> (FLINT_BITS / 4)) % (UWORD(1) << (FLINT_BITS / 4));
            c2 = (c >> (2 * (FLINT_BITS / 4))) % (UWORD(1) << (FLINT_BITS / 4));
            c3 = (c >> (3 * (FLINT_BITS / 4)));

            res[i + 0] = nmod_set_ui(c0, mod);
            res[i + 1] = nmod_set_ui(c1, mod);
            res[i + 2] = nmod_set_ui(c2, mod);
            res[i + 3] = nmod_set_ui(c3, mod);
        }

        if (n % 4)
        {
            c = 0;
            for (j = 0; j < blen; j++)
                c += apre[(i / 4) * n + j] * b[j];

            c0 = c % (UWORD(1) << (FLINT_BITS / 4));
            res[i + 0] = nmod_set_ui(c0, mod);

            if (n % 4 >= 2)
            {
                c1 = (c >> (FLINT_BITS / 4)) % (UWORD(1) << (FLINT_BITS / 4));
                res[i + 1] = nmod_set_ui(c1, mod);

                if (n % 4 == 3)
                {
                    c2 = (c >> (2 * (FLINT_BITS / 4)));
                    res[i + 2] = nmod_set_ui(c2, mod);
                }
            }
        }
    }
}

void
_nmod_poly_mulmod_precond_init_method(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, int method, nmod_t mod)
{
    if (method == NMOD_POLY_MULMOD_PRECOND_NONE)
        _nmod_poly_mulmod_precond_init_none(precond, a, alen, d, dlen, dinv, lendinv, mod);
    else if (method == NMOD_POLY_MULMOD_PRECOND_MATRIX)
        _nmod_poly_mulmod_precond_init_matrix(precond, a, alen, d, dlen, dinv[0], mod);
    else if (method == NMOD_POLY_MULMOD_PRECOND_SHOUP)
        _nmod_poly_mulmod_precond_init_shoup(precond, a, alen, d, dlen, dinv, lendinv, mod);
    else
        flint_throw(FLINT_ERROR, "_nmod_poly_mulmod_precond_init_method: invalid method\n");
}

void
nmod_poly_mulmod_precond_init_method(nmod_poly_mulmod_precond_t precond, const nmod_poly_t a, const nmod_poly_t d, const nmod_poly_t dinv, int method)
{
    _nmod_poly_mulmod_precond_init_method(precond, a->coeffs, a->length, d->coeffs, d->length, dinv->coeffs, dinv->length, method, d->mod);
}

void
_nmod_poly_mulmod_precond_init_num(nmod_poly_mulmod_precond_t precond, nn_srcptr a, slong alen, nn_srcptr d, slong dlen, nn_srcptr dinv, slong lendinv, slong num, nmod_t mod)
{
    int method = _nmod_poly_mulmod_precond_select_method(d, dlen, num, mod);

    _nmod_poly_mulmod_precond_init_method(precond, a, alen, d, dlen, dinv, lendinv, method, mod);
}

void
nmod_poly_mulmod_precond_init_num(nmod_poly_mulmod_precond_t precond, const nmod_poly_t a, const nmod_poly_t d, const nmod_poly_t dinv, slong num)
{
    _nmod_poly_mulmod_precond_init_num(precond, a->coeffs, a->length, d->coeffs, d->length, dinv->coeffs, dinv->length, num, d->mod);
}

void
nmod_poly_mulmod_precond_clear(nmod_poly_mulmod_precond_t precond)
{
    if (precond->method == NMOD_POLY_MULMOD_PRECOND_SHOUP)
        flint_free(precond->adivd);
    if (precond->method == NMOD_POLY_MULMOD_PRECOND_MATRIX)
        flint_free(precond->matrix);
}


void _nmod_poly_mulmod_precond(nn_ptr res, const nmod_poly_mulmod_precond_t apre, nn_srcptr b, slong blen, nmod_t mod)
{
    if (apre->method == NMOD_POLY_MULMOD_PRECOND_NONE)
    {
        slong alen = apre->alen;

        if (alen + blen - 1 <= apre->n)
        {
            if (alen >= blen)
                _nmod_poly_mul(res, apre->a, alen, b, blen, mod);
            else
                _nmod_poly_mul(res, b, blen, apre->a, alen, mod);
            _nmod_vec_zero(res + alen + blen - 1, apre->n - (alen + blen - 1));
        }
        else
        {
            _nmod_poly_mulmod_preinv(res, apre->a, alen, b, blen, apre->d, apre->n + 1, apre->dinv, apre->lendinv, mod);
        }
    }
    else if (apre->method == NMOD_POLY_MULMOD_PRECOND_MATRIX)
        _nmod_poly_mulmod_precond_matrix(res, apre->matrix, apre->n, b, blen, apre->packing, apre->dot_params, mod);
    else
        _nmod_poly_mulmod_precond_shoup(res, apre->a, apre->alen, apre->adivd, apre->n, b, blen, apre->d, mod);
}

void
nmod_poly_mulmod_precond(nmod_poly_t res, const nmod_poly_mulmod_precond_t apre, const nmod_poly_t b)
{
    slong n = apre->n;

    if (apre->alen == 0 || b->length == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    FLINT_ASSERT(b->length <= n);
    FLINT_ASSERT(res->coeffs != apre->a);

    if (res == b)
    {
        nmod_poly_t tmp;
        nmod_poly_init_mod(tmp, b->mod);
        nmod_poly_fit_length(tmp, n);
        _nmod_poly_mulmod_precond(tmp->coeffs, apre, b->coeffs, b->length, res->mod);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
    }
    else
    {
        nmod_poly_fit_length(res, n);
        _nmod_poly_mulmod_precond(res->coeffs, apre, b->coeffs, b->length, res->mod);
    }

    _nmod_poly_set_length(res, n);
    _nmod_poly_normalise(res);
}

