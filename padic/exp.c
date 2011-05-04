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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "padic.h"

/*
    Computes
    \begin{equation*}
    \exp(x) = \sum_{i=0}^{\infty} \frac{x^i}{i!}
    \end{equation*}
    modulo $p^N$.

    Assumes that $\ord_p(x) > 1 / (p - 1)$ for otherwise the 
    function does not converge.  Let $d = \ord_p(x) - 1 / (p - 1) > 0$.

    Note $\ord_p(i!) \leq (i - 1) / (p - 1)$ and hence 
    \begin{align*}
    \ord_p(x^i / i!) & = i \ord_p(x) - \ord_p(i!) \\
                     & \geq i \ord_p(x) - (i - 1) / (p - 1) \\
                     & = i (\ord_p(x) - 1 / (p - 1)) + 1 / (p - 1) \\
                     & = i d + 1 / (p - 1)
    \end{align*}

    We want $i d + 1 / (p - 1) \geq N$ so $i \geq (N - 1 / (p - 1)) / d$.

    Does not support aliasing between \code{rop} and \code{p}.

    Assumes that $p$ is an odd prime.
 */
void _padic_exp_p(fmpz_t rop, const fmpz_t p, long N)
{
    fmpz_t i, hi;

    if (N <= 2)
    {
        if (N == 1) fmpz_set_ui(rop, 1);
        if (N == 2) fmpz_add_ui(rop, p, 1);
        return;
    }

    fmpz_init(i);
    fmpz_init(hi);

    /* Compute hi */
    {
        fmpz_t num, den;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_mul_si(num, p, N + 1);
        fmpz_sub_ui(num, num, N + 4);
        fmpz_sub_ui(den, p, 2);
        fmpz_fdiv_q(hi, num, den);
        fmpz_clear(num);
        fmpz_clear(den);
    }

    fmpz_add_ui(rop, p, 1);

    {
        fmpz_t u;
        long v;

        fmpz_t j, iinv;

        fmpz_t pow;

        fmpz_init(u);
        fmpz_init(j);
        fmpz_init(iinv);
        fmpz_init(pow);

        fmpz_set_ui(u, 1);
        v = 1;

        for (fmpz_set_ui(i, 2); fmpz_cmp(i, hi) < 0; fmpz_add_ui(i, i, 1))
        {
            v += 1 - fmpz_remove(j, i, p);

            if (v >= N)
                break;

            _padic_inv(iinv, j, p, N - v);

            fmpz_pow_ui(pow, p, N - v);
            fmpz_mul(u, u, iinv);
            fmpz_mod(u, u, pow);

            fmpz_pow_ui(pow, p, v);
            fmpz_addmul(rop, u, pow);
            fmpz_pow_ui(pow, p, N);
            fmpz_mod(rop, rop, pow);
        }

        fmpz_clear(u);
        fmpz_clear(j);
        fmpz_clear(iinv);
        fmpz_clear(pow);
    }

    fmpz_clear(i);
    fmpz_clear(hi);
}

/*
    Returns whether the $p$-adic exponential function converges at 
    the point \code{op}, and if so sets \code{rop} to its value.
 */
int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if (fmpz_cmp_ui(ctx->p, 2) && (op[1] >= 1))
    {
        fmpz_t e, pow;

        fmpz_init(e);
        fmpz_init(pow);

        fmpz_pow_ui(e, ctx->p, op[1] - 1);
        fmpz_mul(e, op, e);
        fmpz_pow_ui(pow, ctx->p, ctx->N);

        _padic_exp_p(rop, ctx->p, ctx->N);

        fmpz_powm(rop, rop, e, pow);

        rop[1] = 0;

        fmpz_clear(e);
        fmpz_clear(pow);

        return 1;
    }
    else if (!(fmpz_cmp_ui(ctx->p, 2)) && (op[1] >= 2))
    {
        printf("Exception (padic_exp).  Case p = 2.\n");
        abort();
    }
    else
    {
        return 0;
    }
}

