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

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "qadic.h"

/* 
    Computes the composition $f(g(X))$ modulo the sparse polynomial 
    given by the data \code{(a, j, lena)}, which is assumed to be 
    of degree~$d \geq 2$.

    Sets the vector \code{(rop, d)}.

    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$.

    Does not support aliasing.
 */

static void 
_fmpz_mod_poly_compose_mod(fmpz *rop, 
                           const fmpz *op1, long len1, 
                           const fmpz *op2, long len2, 
                           const fmpz *a, const long *j, long lena, 
                           const fmpz_t p)
{
    const long d = j[lena - 1];

    if (len1 == 1)
    {
        fmpz_set(rop, op1);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        long i;
        fmpz *t;

        t = _fmpz_vec_init(2*d - 1);

        i = len1 - 1;

        _fmpz_mod_poly_scalar_mul_fmpz(rop, op2, len2, op1 + i, p);
        _fmpz_vec_zero(rop + len2, d - len2);
        i--;
        if (i >= 0)
        {
            fmpz_add(rop, rop, op1 + i);
            if (fmpz_cmpabs(rop, p) >= 0)
                fmpz_sub(rop, rop, p);
        }

        while (i > 0)
        {
            i--;
            _fmpz_mod_poly_mul(t, rop, d, op2, len2, p);
            _fmpz_mod_poly_reduce(t, d + len2 - 1, a, j, lena, p);
            _fmpz_mod_poly_add(rop, t, d, op1 + i, 1, p);
        }

        _fmpz_vec_clear(t, 2*d - 1);
    }
}

/*
    Computes $\sigma(X) \bmod{p^N}$ where $X$ is such that 
    $\mathbf{Q}_q \cong \mathbf{Q}_p[X]/(f(X))$.

    Assumes that the precision $N$ is at least~$2$.

    Sets \code{(rop, 2*d - 1)}.
 */

void _qadic_sigma_a(fmpz *rop, 
                    const fmpz *a, const long *j, long lena, 
                    const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    long *e, i, n;
    fmpz *pow, *f1, *f2, *num, *den, *t;

    n = FLINT_CLOG2(N) + 1;

    e = flint_malloc(n * sizeof(long));
    for (e[i = 0] = N; e[i] > 1; i++)
        e[i + 1] = (e[i] + 1) / 2;

    pow = _fmpz_vec_init(n);
    num = _fmpz_vec_init(d);
    den = _fmpz_vec_init(d);
    f1  = _fmpz_vec_init(d + 1);
    f2  = _fmpz_vec_init(d);
    t   = _fmpz_vec_init(2*d - 1);

    /* Compute powers of p */
    {
        fmpz_one(t);
        fmpz_set(pow + i, p);
    }
    for (i--; i >= 1; i--)
    {
        if (e[i] & 1L)
        {
            fmpz_mul(pow + i, t, pow + (i + 1));
            fmpz_mul(t, t, t);
        }
        else
        {
            fmpz_mul(t, t, pow + (i + 1));
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }
    }
    {
        if (e[i] & 1L)
            fmpz_mul(pow + i, t, pow + (i + 1));
        else
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
    }

    /* Dense representation of f and f' */
    {
        long k;

        for (k = 0; k < lena; k++)
            fmpz_set(f1 + j[k], a + k);
        for (k = 0; k < lena; k++)
            if (j[k] >= 1)
                fmpz_mul_ui(f2 + (j[k] - 1), a + k, j[k]);
    }

    /* Run Newton iteration */
    i = n - 1;
    {
        fmpz op[2] = {0L, 1L};

        _qadic_pow(rop, op, 2, p, a, j, lena, pow + i);
    }
    for (i--; i >= 0; i--)
    {
        /* z' := z - f(z) / f'(z) */

        _fmpz_mod_poly_compose_mod(num, f1, d + 1, rop, d, a, j, lena, pow + i);
        _fmpz_mod_poly_compose_mod(den, f2, d, rop, d, a, j, lena, pow + i);

        _qadic_inv(den, den, d, a, j, lena, p, e[i]);

        _fmpz_mod_poly_mul(t, num, d, den, d, pow + i);
        _fmpz_mod_poly_reduce(t, 2*d - 1, a,j, lena, pow + i);
        _fmpz_mod_poly_sub(rop, rop, d, t, d, pow + i);
    }

    _fmpz_vec_clear(pow, n);
    _fmpz_vec_clear(num, d);
    _fmpz_vec_clear(den, d);
    _fmpz_vec_clear(f1, d + 1);
    _fmpz_vec_clear(f2, d);
    _fmpz_vec_clear(t, 2*d - 1);
    flint_free(e);
}

/*
    Sets \code{(rop, d)} but requires \code{rop} to be an array of 
    length at least $2d - 1$.
 */

void _qadic_sigma(fmpz *rop, const fmpz *op, long len, 
                  const fmpz *a, const long *j, long lena, 
                  const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    if (len == 1)  /* op is in Zp, not just Zq */
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, d - len);
    }
    else if (N == 1)
    {
        _qadic_pow(rop, op, len, p, a, j, lena, p);
    }
    else
    {
        fmpz *t;
        fmpz_t pow;

        t = _fmpz_vec_init(2*d - 1);
        fmpz_init(pow);
        fmpz_pow_ui(pow, p, N);

        _qadic_sigma_a(t, a, j, lena, p, N);

        _fmpz_mod_poly_compose_mod(rop, op, len, t, d, a, j, lena, pow);
        _fmpz_vec_zero(rop + d, d - 1);

        _fmpz_vec_clear(t, 2*d - 1);
        fmpz_clear(pow);
    }
}

void qadic_sigma(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (qadic_is_zero(op) || op->val >= N)
    {
        qadic_zero(rop);
    }
    else
    {
        const long d = qadic_ctx_degree(ctx);
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            padic_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _qadic_sigma(t, op->coeffs, op->length, 
                     ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N - op->val);

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            rop->val = op->val;
            _padic_poly_set_length(rop, d);
        }
        _padic_poly_normalise(rop);
    }
}

