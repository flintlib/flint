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
    Copyright (C) 2012 Fredrik Johansson
 
******************************************************************************/

#include "padic.h"
#include "ulong_extras.h"

static void
_padic_log_bsplit_series(fmpz_t P, fmpz_t B, fmpz_t T, 
                         const fmpz_t x, len_t a, len_t b)
{
    if (b - a == 1)
    {
        fmpz_set(P, x);
        fmpz_set_si(B, a);
        fmpz_set(T, x);
    }
    else if (b - a == 2)
    {
        fmpz_mul(P, x, x);
        fmpz_set_si(B, a);
        fmpz_mul_si(B, B, a + 1);
        fmpz_mul_si(T, x, a + 1);
        fmpz_addmul_ui(T, P, a);
    }
    else
    {
        const len_t m = (a + b) / 2;

        fmpz_t RP, RB, RT;

        _padic_log_bsplit_series(P, B, T, x, a, m);

        fmpz_init(RP);
        fmpz_init(RB);
        fmpz_init(RT);

        _padic_log_bsplit_series(RP, RB, RT, x, m, b);

        fmpz_mul(RT, RT, P);
        fmpz_mul(T, T, RB);
        fmpz_addmul(T, RT, B);
        fmpz_mul(P, P, RP);
        fmpz_mul(B, B, RB);

        fmpz_clear(RP);
        fmpz_clear(RB);
        fmpz_clear(RT);
    }
}

/*
    Assumes that $y = 1 - x$ is such that $\log(x)$ 
    converges.

    Assumes that $v = \ord_p(y)$ with $v < N$, which 
    also forces $N$ to be positive.

    The result $z$ might not be reduced modulo $p^N$.

    Supports aliasing between $y$ and $z$.
 */

static void 
_padic_log_bsplit(fmpz_t z, const fmpz_t y, len_t v, const fmpz_t p, len_t N)
{
    fmpz_t P, B, T;
    len_t k, n;

    n = _padic_log_bound(v, N, p);
    n = FLINT_MAX(n, 2);

    fmpz_init(P);
    fmpz_init(B);
    fmpz_init(T);

    _padic_log_bsplit_series(P, B, T, y, 1, n);

    k = fmpz_remove(B, B, p);
    fmpz_pow_ui(P, p, k);
    fmpz_divexact(T, T, P);

    _padic_inv(B, B, p, N);
    fmpz_mul(z, T, B);

    fmpz_clear(P);
    fmpz_clear(B);
    fmpz_clear(T);
}

void 
_padic_log_balanced(fmpz_t z, const fmpz_t y, len_t v, const fmpz_t p, len_t N)
{
    fmpz_t pv, pN, r, t, u;
    len_t w;
    padic_inv_t S;

    fmpz_init(pv);
    fmpz_init(pN);
    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(u);
    _padic_inv_precompute(S, p, N);

    fmpz_set(pv, p);
    fmpz_pow_ui(pN, p, N);
    fmpz_mod(t, y, pN);
    fmpz_zero(z);
    w = 1;

    while (!fmpz_is_zero(t))
    {
        fmpz_mul(pv, pv, pv);
        fmpz_fdiv_qr(t, r, t, pv);

        if (!fmpz_is_zero(t))
        {
            fmpz_mul(t, t, pv);
            fmpz_sub_ui(u, r, 1);
            fmpz_neg(u, u);
            _padic_inv_precomp(u, u, S);
            fmpz_mul(t, t, u);
            fmpz_mod(t, t, pN);
        }

        if (!fmpz_is_zero(r))
        {
            _padic_log_bsplit(r, r, w, p, N);
            fmpz_sub(z, z, r);
        }
        w *= 2;
    }

    fmpz_mod(z, z, pN);

    fmpz_clear(pv);
    fmpz_clear(pN);
    fmpz_clear(r);
    fmpz_clear(t);
    fmpz_clear(u);
    _padic_inv_clear(S);
}

int padic_log_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const fmpz *p = ctx->p;
    const len_t N  = padic_prec(rop);

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
            len_t v;

            fmpz_init(t);
            v = fmpz_remove(t, x, p);
            fmpz_clear(t);

            if (v >= 2 || (!fmpz_equal_ui(p, 2) && v >= 1))
            {
                if (v >= N)
                {
                    padic_zero(rop);
                }
                else
                {
                    _padic_log_balanced(padic_unit(rop), x, v, p, N);
                    padic_val(rop) = 0;
                    _padic_canonicalise(rop, ctx);
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

