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
 */

/*
    Assumes that $p$ is an odd prime.

    The above assumptions simplify and defining 
    $M = N + \floor{(N - 2)/(p - 2)}$ imply that 
    need to consider $i$ in the range $[0, M]$.

    Does not support aliasing between \code{rop} and \code{p}.
 */
void _padic_exp_p(fmpz_t rop, const fmpz_t p, long N)
{
    fmpz_t t, f, pow;
    fmpz *j;
    long i, M, n;

    if (N <= 2)
    {
        if (N == 1) fmpz_set_ui(rop, 1);
        if (N == 2) fmpz_add_ui(rop, p, 1);
        return;
    }

    fmpz_init(t);
    fmpz_init(f);
    fmpz_init(pow);

    M = N + (fmpz_fits_si(p) ? (N - 2) / (fmpz_get_si(p) - 2) : 0);

    fmpz_set_si(t, M);
    padic_val_fac(t, t, p);
    n = fmpz_get_si(t);

    fmpz_pow_ui(pow, p, N + n);

    /* j[M] = 1, j[M-i] = M ... (M-(i+1)) */
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

    fmpz_clear(t);
    fmpz_clear(f);
    fmpz_clear(pow);
    _fmpz_vec_clear(j, M + 1);
}

/*
    Evalutes $\exp_2(4) \pmod{2^N}$.

    This simplifies as follows:

    Let $h = N - 2$.  Then we compute 
    \begin{equation*}
    \frac{1}{h!} \sum_{i=0}^{h} 4^i (i + 1) \dotsm h \pmod{2^N}
    \end{equation*}
    as 
    \begin{equation*}
    \frac{2^{\ord_2(h!)}}{h!} \sum_{i=0}^{h} 2^{2i - \ord_2(h!)} (i + 1) \dotsm h \pmod{2^N}
    \end{equation*}
 */
void _padic_exp_4(fmpz_t rop, long N)
{
    if (N <= 2)
    {
        fmpz_set_ui(rop, 1);
    }
    else
    {
        long i, n, h = N - 2;
        fmpz_t p, s, t, u;

        fmpz_init(p);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(u);

        n = padic_val_fac_ui2(h);

        fmpz_set_ui(p, 2);
        fmpz_set_ui(s, 1);
        fmpz_set_ui(t, 1);
        fmpz_mul_2exp(u, t, 2 * h - n);
        fmpz_fdiv_r_2exp(rop, u, N);

        for (i = h; i > 0; i--)
        {
            int c;

            count_trailing_zeros(c, i);

            fmpz_mul_ui(s, s, i >> c);
            fmpz_fdiv_r_2exp(s, s, N);

            fmpz_mul_ui(t, t, i);
            fmpz_fdiv_r_2exp(t, t, N + n);

            if (2 * (i - 1) - n > 0)
                fmpz_mul_2exp(u, t, 2 * (i - 1) - n);
            else
                fmpz_fdiv_q_2exp(u, t, - 2 * (i - 1) + n);
            fmpz_fdiv_r_2exp(u, u, N);

            fmpz_add(rop, rop, u);
            if (i % 4)
                fmpz_fdiv_r_2exp(rop, rop, N);
        }

        _padic_inv(s, s, p, N);

        fmpz_mul(rop, rop, s);
        fmpz_fdiv_r_2exp(rop, rop, N);

        fmpz_clear(p);
        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(u);
    }
}

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if ((fmpz_cmp_ui(ctx->p, 2) == 0) && (padic_val(op) >= 2))
    {
        fmpz_t e, m;

        fmpz_init(e);
        fmpz_init(m);

        fmpz_pow_ui(e, ctx->p, padic_val(op) - 2);
        fmpz_mul(e, padic_unit(op), e);

        _padic_exp_4(padic_unit(rop), ctx->N);

        fmpz_set_ui(m, 1);
        fmpz_mul_2exp(m, m, ctx->N);
        fmpz_powm(padic_unit(rop), padic_unit(rop), e, m);
        padic_val(rop) = 0;

        fmpz_clear(m);
        fmpz_clear(e);

        return 1;
    }
    else if ((fmpz_cmp_ui(ctx->p, 2) > 0) && (padic_val(op) >= 1))
    {
        int alloc;
        fmpz_t e, m;

        fmpz_init(e);
        fmpz_pow_ui(e, ctx->p, padic_val(op) - 1);
        fmpz_mul(e, padic_unit(op), e);

        _padic_exp_p(padic_unit(rop), ctx->p, ctx->N);

        padic_val(rop) = 0;
        _padic_ctx_pow_ui(m, &alloc, ctx->N, ctx);
        fmpz_powm(padic_unit(rop), padic_unit(rop), e, m);
        if (alloc)
            fmpz_clear(m);
        fmpz_clear(e);

        return 1;
    }
    else
    {
        return 0;
    }
}

