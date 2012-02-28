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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"

/*
    Returns $b$ such that for all $i \geq b$ we have 
    \begin{equation*}
    i v - \ord_p(i) \geq N
    \end{equation*}
    where $v \geq 1$.

    Assumes that $1 \leq v < N$.

    Assumes that $b v$, and $N + e$ do not overflow, 
    where $e = \floor{\log_{p}{b}}$.
 */
long _padic_log_bound(long v, long N, long p)
{
    long e, i = (N - 1) / v;
    mp_limb_t j;

    do 
    {
        j = ++i;
        e = n_remove(&j, p);
    }
    while (i * v < N + e);

    return i;
}

/*
    Computes 
    \begin{equation*}
    z = \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    \end{equation*}

    Note that this can be used to compute the $p$-adic logarithm 
    via the equation 
    \begin{align*}
    \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
            & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    \end{align*}

    Assumes that $y = 1 - x$ is non-zero and that $v = \ord_p(y)$ 
    is at least $1$ when $p$ is odd and at least $2$ when $p = 2$ 
    so that the series converges.

    Assumes that $v < N$.

    Does not support aliasing between $y$ and $z$.
 */
static void _padic_log(fmpz_t z, const fmpz_t y, long v, const padic_ctx_t ctx)
{
    if (fmpz_fits_si(ctx->p))
    {
        /*
            Assumes that the index i fits into a small fmpz.
         */
        long e, i, j, k, p;
        fmpz_t m, s, t;
        fmpz *q;
        padic_inv_t pre;

        p = fmpz_get_si(ctx->p);
        i = _padic_log_bound(v, ctx->N, p) - 1;

        k = n_flog(i, p);

        fmpz_init(m);
        fmpz_init(s);
        fmpz_init(t);
        q = _fmpz_vec_init(k + 1);

        _padic_inv_precompute(pre, ctx->p, ctx->N + k);

        fmpz_pow_ui(m, ctx->p, ctx->N + k);
        fmpz_one(q + 0);
        for (j = 1; j <= k; j++)
            fmpz_mul_ui(q + j, q + (j - 1), p);

        fmpz_zero(z);

        for ( ; i > 0; i--)
        {
            fmpz_mul(t, z, y);

            j = i;
            e = n_remove((mp_limb_t *) &j, p);
            _padic_inv_precomp(s, (fmpz *) &j, pre);
            fmpz_mul(z, s, q + (k - e));

            fmpz_add(z, z, t);
            fmpz_mod(z, z, m);
        }

        fmpz_divexact(z, z, q + k);
        fmpz_mul(z, z, y);

        fmpz_clear(m);
        fmpz_clear(s);
        fmpz_clear(t);
        _fmpz_vec_clear(q, k + 1);
        _padic_inv_clear(pre);
    }
    else
    {
        /*
            When p does not fit into a signed long, 
            p does not divide the index i.

            Assumes that (N - 1) / v is a small 
            fmpz integer.
         */
        long i;
        fmpz_t m, t;

        i = (ctx->N - 1) / v;

        fmpz_init(m);
        fmpz_init(t);

        fmpz_pow_ui(m, ctx->p, ctx->N);

        fmpz_zero(z);

        for ( ; i > 0; i--)
        {
            fmpz_mul(t, z, y);

            _padic_inv(z, (fmpz *) &i, ctx->p, ctx->N);

            fmpz_add(z, z, t);
            fmpz_mod(z, z, m);
        }

        fmpz_mul(z, z, y);

        fmpz_clear(m);
        fmpz_clear(t);
    }
}

int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_val(op) < 0)
    {
        return 0;
    }
    else
    {
        fmpz_t x;
        int ans;

        fmpz_init(x);

        padic_get_fmpz(x, op, ctx);
        fmpz_sub_ui(x, x, 1);
        fmpz_neg(x, x);

        if (fmpz_is_zero(x))
        {
            padic_zero(rop);
            ans = 1;
        }
        else
        {
            fmpz_t t;
            long v;

            fmpz_init(t);
            v = fmpz_remove(t, x, ctx->p);
            fmpz_clear(t);

            if ((*(ctx->p) == 2L && v >= 2) || v >= 1)
            {
                if (v >= ctx->N)
                {
                    padic_zero(rop);
                }
                else
                {
                    _padic_log(padic_unit(rop), x, v, ctx);
                    fmpz_neg(padic_unit(rop), padic_unit(rop));
                    padic_val(rop) = 0;
                    padic_reduce(rop, ctx);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        fmpz_clear(x);
        return ans;
    }
}

