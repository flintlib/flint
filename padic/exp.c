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

    Assumes that $\ord_p(x) > 1 / (p - 1)$ for otherwise the function 
    does not converge.  Let $d = \ord_p(x) - 1 / (p - 1) > 0$.

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
    fmpz_t hi, t, f, pow;
    fmpz *j;
    long i, M, n;

    if (N <= 2)
    {
        if (N == 1) fmpz_set_ui(rop, 1);
        if (N == 2) fmpz_add_ui(rop, p, 1);
        return;
    }

    fmpz_init(hi);
    fmpz_init(t);
    fmpz_init(f);
    fmpz_init(pow);

    /* Compute hi */
    {
        fmpz_t num, den;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_mul_si(num, p, N + 1);
        fmpz_sub_ui(num, num, N + 4);
        fmpz_sub_ui(den, p, 2);
        fmpz_fdiv_q(hi, num, den);
        fmpz_sub_ui(hi, hi, 1);
        fmpz_clear(num);
        fmpz_clear(den);
    }

    if (COEFF_IS_MPZ(*hi))
    {
        printf("Exception (_padic_exp_p).\n");
        printf("Computing the p-adic exponential requires too many terms.\n");
        abort();
    }

    M = *hi;

    padic_val_fac(t, hi, p);
    n = *t;

    fmpz_pow_ui(pow, p, N + n);

    j = _fmpz_vec_init(M + 1);

    fmpz_set_ui(j + M, 1);
    for (i = M; i > 0; i--)
    {
        fmpz_mul_ui(j + (i - 1), j + i, i);
        fmpz_mod(j + (i - 1), j + (i - 1), pow);
    }

    fmpz_add_ui(rop, p, 1);
    fmpz_mul(rop, rop, j + 0);

    fmpz_set(f, p);
    for (i = 2; i <= M; i++)
    {
        fmpz_mul(f, f, p);
        fmpz_mul(t, f, j + i);
        fmpz_add(rop, rop, t);
        if (i & 1L)
            fmpz_mod(rop, rop, pow);
    }

    _padic_inv(t, j + 0, p, N + n);

    fmpz_mul(rop, rop, t);
    fmpz_mod(rop, rop, pow);

    fmpz_clear(hi);
    fmpz_clear(t);
    fmpz_clear(f);
    fmpz_clear(pow);
    _fmpz_vec_clear(j, M + 1);
}

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if (fmpz_cmp_ui(ctx->p, 2) && (padic_val(op)>= 1))
    {
        fmpz_t e, m, pow;
        int alloc;

        fmpz_init(e);
        fmpz_init(m);

        fmpz_pow_ui(e, ctx->p, padic_val(op) - 1);
        fmpz_mul(e, padic_unit(op), e);
        fmpz_sub_ui(m, ctx->p, 1);
        _padic_ctx_pow_ui(pow, &alloc, ctx->N - 1, ctx);
        fmpz_mul(m, m, pow);
        if (alloc) 
            fmpz_clear(pow);
        fmpz_mod(e, e, m);

        _padic_exp_p(padic_unit(rop), ctx->p, ctx->N);

        _padic_ctx_pow_ui(pow, &alloc, ctx->N, ctx);
        fmpz_powm(padic_unit(rop), padic_unit(rop), e, pow);
        if (alloc)
            fmpz_clear(pow);

        padic_val(rop) = 0;

        fmpz_clear(e);
        fmpz_clear(m);

        return 1;
    }
    else if (!(fmpz_cmp_ui(ctx->p, 2)) && (padic_val(op) >= 2))
    {
        printf("Exception (padic_exp).  Case p = 2 not implemented yet.\n");
        abort();
    }
    else
    {
        return 0;
    }
}

