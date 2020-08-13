/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qadic.h"

/*
    Computes the characteristic polynomial of the $n \times n$ matrix $M$ 
    modulo \code{pN} using a division-free algorithm in $O(n^4)$ ring 
    operations.

    Only returns the determinant.

    Assumes that $n$ is at least $2$.
 */

static 
void _fmpz_mod_mat_det(fmpz_t rop, const fmpz *M, slong n, const fmpz_t pN)
{
    fmpz *F;
    fmpz *a;
    fmpz *A;
    fmpz_t s;
    slong t, i, j, p, k;

    F = _fmpz_vec_init(n);
    a = _fmpz_vec_init((n-1) * n);
    A = _fmpz_vec_init(n);

    fmpz_init(s);

    fmpz_neg(F + 0, M + 0*n + 0);

    for (t = 1; t < n; t++)
    {
        for (i = 0; i <= t; i++)
            fmpz_set(a + 0*n + i, M + i*n + t);

        fmpz_set(A + 0, M + t*n + t);

        for (p = 1; p < t; p++)
        {
            for (i = 0; i <= t; i++)
            {
                fmpz_zero(s);
                for (j = 0; j <= t; j++)
                    fmpz_addmul(s, M + i*n + j, a + (p-1)*n + j);
                fmpz_mod(a + p*n + i, s, pN);
            }

            fmpz_set(A + p, a + p*n + t);
        }

        fmpz_zero(s);
        for (j = 0; j <= t; j++)
            fmpz_addmul(s, M + t*n + j, a + (t-1)*n + j);
        fmpz_mod(A + t, s, pN);

        for (p = 0; p <= t; p++)
        {
            fmpz_sub(F + p, F + p, A + p);
            for (k = 0; k < p; k++)
                fmpz_submul(F + p, A + k, F + (p-k-1));
            fmpz_mod(F + p, F + p, pN);
        }
    }

    /*
        Now [F{n-1}, F{n-2}, ..., F{0}, 1] is the 
        characteristic polynomial of the matrix M.
     */

    if (n % WORD(2) == 0)
    {
        fmpz_set(rop, F + (n-1));
    }
    else
    {
        fmpz_neg(rop, F + (n-1));
        fmpz_mod(rop, rop, pN);
    }

    _fmpz_vec_clear(F, n);
    _fmpz_vec_clear(a, (n-1)*n);
    _fmpz_vec_clear(A, n);
    fmpz_clear(s);
}

void _qadic_norm_resultant(fmpz_t rop, const fmpz *op, slong len, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    fmpz_t pN;

    fmpz_init(pN);
    fmpz_pow_ui(pN, p, N);

    if (len == 1)
    {
        fmpz_powm_ui(rop, op + 0, d, pN);
    }
    else  /* len >= 2 */
    {
        {
            const slong n = d + len - 1;
            slong i, k;
            fmpz *M;

            M = flint_calloc(n * n, sizeof(fmpz));

            for (k = 0; k < len-1; k++)
            {
                for (i = 0; i < lena; i++)
                {
                    M[k * n + k + (d - j[i])] = a[i];
                }
            }
            for (k = 0; k < d; k++)
            {
                for (i = 0; i < len; i++)
                {
                    M[(len-1 + k) * n + k + (len-1 - i)] = op[i];
                }
            }

            _fmpz_mod_mat_det(rop, M, n, pN);
            
            flint_free(M);
        }

        /*
            XXX:  This part of the code is currently untested as the Conway 
            polynomials used for the extension Qq/Qp are monic.
         */
        if (!fmpz_is_one(a + (lena - 1)))
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_powm_ui(f, a + (lena - 1), len - 1, pN);
            _padic_inv(f, f, p, N);
            fmpz_mul(rop, f, rop);
            fmpz_mod(rop, rop, pN);
            fmpz_clear(f);
        }
    }
    fmpz_clear(pN);
}

void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const slong N = padic_prec(rop);
    const slong d = qadic_ctx_degree(ctx);

    /* N(p^v u) = p^{dv} N(u) */

    if (qadic_is_zero(op) || d * op->val >= N)
    {
        padic_zero(rop);
    }
    else
    {
        _qadic_norm_resultant(padic_unit(rop), op->coeffs, op->length, 
                              ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, 
                              N - d * op->val);
        padic_val(rop) = d * op->val;
    }
}

