/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

void
_nmod_mat_charpoly_berkowitz(nn_ptr cp, const nmod_mat_t mat, nmod_t mod)
{
    const slong n = mat->r;

    if (mod.n == 1)
    {
        _nmod_vec_zero(cp, n + 1);
    }
    else if (n == 0)
    {
        cp[0] = 1;
    }
    else if (n == 1)
    {
        cp[0] = nmod_neg(nmod_mat_entry(mat, 0, 0), mod);
        cp[1] = 1;
    }
    else if (n == 2)
    {
        cp[0] = nmod_sub(nmod_mul(nmod_mat_entry(mat, 0, 0), nmod_mat_entry(mat, 1, 1), mod),
                         nmod_mul(nmod_mat_entry(mat, 0, 1), nmod_mat_entry(mat, 1, 0), mod), mod);
        cp[1] = nmod_add(nmod_mat_entry(mat, 0, 0), nmod_mat_entry(mat, 1, 1), mod);
        cp[1] = nmod_neg(cp[1], mod);
        cp[2] = 1;
    }
    else
    {
        slong i, k, t;
        nn_ptr a, A, s;
        TMP_INIT;

        TMP_START;
        a = TMP_ALLOC(sizeof(ulong) * (n * n));
        A = a + (n - 1) * n;

        const dot_params_t params = _nmod_vec_dot_params(n, mod);

        _nmod_vec_zero(cp, n + 1);
        cp[0] = nmod_neg(nmod_mat_entry(mat, 0, 0), mod);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                a[0 * n + i] = nmod_mat_entry(mat, i, t);
            }

            A[0] = nmod_mat_entry(mat, t, t);

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    s[0] = _nmod_vec_dot(nmod_mat_entry_ptr(mat, i, 0), a + (k - 1) * n, t + 1, mod, params);
                }

                A[k] = a[k * n + t];
            }

            A[t] = _nmod_vec_dot(nmod_mat_entry_ptr(mat, t, 0), a + (t - 1) * n, t + 1, mod, params);

            for (k = 0; k <= t; k++)
            {
                cp[k] = nmod_sub(cp[k], _nmod_vec_dot_rev(A, cp, k, mod, params), mod);
                cp[k] = nmod_sub(cp[k], A[k], mod);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            cp[i] = cp[i - 1];

        cp[0] = 1;
        _nmod_poly_reverse(cp, cp, n + 1, n + 1);

        TMP_END;
    }
}

void nmod_mat_charpoly_berkowitz(nmod_poly_t cp, const nmod_mat_t mat)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly_berkowitz).  Non-square matrix.\n");
    }

    nmod_poly_fit_length(cp, mat->r + 1);
    _nmod_poly_set_length(cp, mat->r + 1);
    _nmod_mat_charpoly_berkowitz(cp->coeffs, mat, mat->mod);
}

static
int _nmod_mat_charpoly_danilevsky(nn_ptr p, const nmod_mat_t M, nmod_t mod)
{
    slong n = M->r, i, j, k, plen;
    ulong * V, * W, * b, * t;
    ulong g, h;
    nmod_mat_t M2;
    int success = 1;
    TMP_INIT;


    if (n == 0)
    {
        p[0] = 1;
        return 1;
    }

    if (n == 1)
    {
        p[1] = 1;
        p[0] = nmod_neg(nmod_mat_entry(M, 0, 0), mod);
        return 1;
    }

    TMP_START;

    i = 1;
    plen = 1;
    p[0] = 1;

    const dot_params_t params = _nmod_vec_dot_params(n, mod);

    nmod_mat_init_set(M2, M);
    b = (ulong *) TMP_ALLOC((n + 1)*sizeof(ulong));
    t = (ulong *) TMP_ALLOC((n + 1)*sizeof(ulong));
    V = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    W = (ulong *) TMP_ALLOC(n*sizeof(ulong));

#define A(ii, jj) nmod_mat_entry(M2, ii, jj)
#define Aptr(ii, jj) nmod_mat_entry_ptr(M2, ii, jj)

    while (i < n)
    {
        h = A(n - i, n - i - 1);

        while (h == 0)
        {
            k = 1;
            while (k < n - i && A(n - i, n - i - k - 1) == 0)
                k++;

            if (k == n - i)
            {
                b[i] = 1;
                for (k = 1; k <= i; k++)
                    b[k - 1] = nmod_neg(A(n - i, n - k), mod);
                if (plen >= i + 1)
                    _nmod_poly_mul(t, p, plen, b, i + 1, mod);
                else
                    _nmod_poly_mul(t, b, i + 1, p, plen, mod);
                plen += i;
                _nmod_vec_set(p, t, plen);
                n -= i;
                i = 1;

                if (n == 1)
                {
                    b[1] = 1;
                    b[0] = nmod_neg(A(0, 0), mod);

                    if (plen >= 2)
                        _nmod_poly_mul(t, p, plen, b, 2, mod);
                    else
                        _nmod_poly_mul(t, b, 2, p, plen, mod);

                    plen += 1;
                    _nmod_vec_set(p, t, plen);
                   goto cleanup;
                }
            }
            else
            {
                nmod_mat_swap_rows(M2, NULL, n - i - k - 1, n - i - 1);

                for (j = 1; j <= n - i + 1; j++)
                    FLINT_SWAP(ulong, A(j - 1, n - i - k - 1), A(j - 1, n - i - 1));
            }

            h = A(n - i, n - i - 1);
        }

        // h = nmod_inv(nmod_neg(h, mod), mod);
        h = nmod_neg(h, mod);
        g = n_gcdinv(&h, h, mod.n);
        if (g != 1)
        {
            success = 0;
            goto cleanup;
        }

        _nmod_vec_set(W, Aptr(n - i, 0), n);
        _nmod_vec_scalar_mul_nmod(V, Aptr(n - i, 0), n, h, mod);
        h = nmod_neg(h, mod);

        for (j = 1; j <= n - i; j++)
        {
            ulong c = A(j - 1, n - i - 1);
            nn_ptr row = nmod_mat_row_ptr(M2, j - 1);

            _nmod_vec_scalar_addmul_nmod(row, V, n - i - 1, c, mod);
            _nmod_vec_scalar_addmul_nmod(row + (n - i), V + (n - i), i, c, mod);

            A(j - 1, n - i - 1) = nmod_mul(c, h, mod);
        }

        for (j = 1; j <= n - i - 1; j++)
            A(n - i - 1, j - 1) = _nmod_vec_dot_strided(Aptr(0, j - 1), M2->stride, W, 1, n - i, mod, params);

        for (j = n - i; j <= n - 1; j++)
            A(n - i - 1, j - 1) = nmod_add(_nmod_vec_dot_strided(Aptr(0, j - 1), M2->stride, W, 1, n - i, mod, params), W[j], mod);

        A(n - i - 1, n - 1) = _nmod_vec_dot_strided(Aptr(0, n - 1), M2->stride, W, 1, n - i, mod, params);

        i++;
    }


    b[n] = 1;
    for (i = 1; i <= n; i++)
        b[i - 1] = nmod_neg(A(0, n - i), mod);
    if (plen >= n + 1)
        _nmod_poly_mul(t, p, plen, b, n + 1, mod);
    else
        _nmod_poly_mul(t, b, n + 1, p, plen, mod);
    plen += n;
    _nmod_vec_set(p, t, plen);

#undef A

cleanup:

    nmod_mat_clear(M2);
    TMP_END;

    return success;
}

int nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t mat)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly_danilevsky).  Non-square matrix.\n");
    }

    nmod_poly_fit_length(p, mat->r + 1);
    _nmod_poly_set_length(p, mat->r + 1);
    return _nmod_mat_charpoly_danilevsky(p->coeffs, mat, mat->mod);
}

#include "gr.h"
#include "gr_mat.h"
#include <stdint.h>

void
nmod_mat_charpoly(nmod_poly_t cp, const nmod_mat_t mat)
{
    if (mat->r >= 32 && mat->mod.n <= 255 && n_is_prime(mat->mod.n))
    {
        slong n = mat->r;
        gr_ctx_t ctx;
        gr_ctx_init_nmod8(ctx, mat->mod.n);

        gr_mat_t M;
        uint8_t * t;

        gr_mat_init(M, n, n, ctx);
        t = GR_TMP_ALLOC(n + 1);

        slong i, j;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                ((uint8_t *) (M->entries))[i * M->stride + j] = mat->entries[i * mat->stride + j];

        GR_MUST_SUCCEED(_gr_mat_charpoly(t, M, ctx));

        nmod_poly_fit_length(cp, n + 1);
        _nmod_poly_set_length(cp, n + 1);

        for (i = 0; i <= n; i++)
            cp->coeffs[i] = t[i];

        GR_TMP_FREE(t, n + 1);
        gr_mat_clear(M, ctx);
        return;
    }

    /* Todo: modulus-dependent cutoff */
    if (mat->r > 8 && nmod_mat_charpoly_danilevsky(cp, mat))
        return;

    nmod_mat_charpoly_berkowitz(cp, mat);
}

