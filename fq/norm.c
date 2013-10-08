/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz
 
******************************************************************************/

#include "fq.h"

/*
    Computes the characteristic polynomial of the $n \times n$ matrix $M$ 
    modulo \code{pN} using a division-free algorithm in $O(n^4)$ ring 
    operations.

    Only returns the determinant.

    Assumes that $n$ is at least $2$.
 */

static 
void _fmpz_mod_mat_det(fmpz_t rop, const fmpz *M, long n, const fmpz_t pN)
{
    if (n == 1)
    {
        fmpz_set(rop, M);
    }
    else
    {
        fmpz *F;
        fmpz *a;
        fmpz *A;
        fmpz_t s;
        long t, i, j, p, k;

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

        if (n % 2L == 0)
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
}

/*
    Computes the norm on $\mathbf{Q}_q$ to precision $N \geq 1$. 
    When $N = 1$, this computes the norm on $\mathbf{F}_q$.
 */

void _fq_norm(fmpz_t rop, const fmpz *op, long len, 
                          const fmpz *a, const long *j, long lena, 
                          const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    fmpz *pN;

    if (N == 1)
    {
        pN = (fmpz *) p;  /* XXX:  Read-only */
    }
    else
    {
        pN = flint_malloc(sizeof(fmpz));
        fmpz_init(pN);
        fmpz_pow_ui(pN, p, N);
    }

    if (len == 1)
    {
        fmpz_powm_ui(rop, op + 0, d, pN);
    }
    else 
    {
        /* Construct an ad hoc matrix M and set rop to det(M) */
        {
            const long n = d + len - 1;
            long i, k;
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
            polynomials used for the extension Fq/Fp are monic.
         */
        if (!fmpz_is_one(a + (lena - 1)))
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_powm_ui(f, a + (lena - 1), len - 1, pN);
            fmpz_invmod(f, f, pN);
            fmpz_mul(rop, f, rop);
            fmpz_mod(rop, rop, pN);
            fmpz_clear(f);
        }
    }

    if (N > 1)
    {
        fmpz_clear(pN);
        flint_free(pN);
    }
}

void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
{
    if (fq_is_zero(op))
    {
        fmpz_zero(rop);
    }
    else
    {
        _fq_norm(rop, op->coeffs, op->length, 
                 ctx->a, ctx->j, ctx->len, fq_ctx_prime(ctx), 1);
    }
}

