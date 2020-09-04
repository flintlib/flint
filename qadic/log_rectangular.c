/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "qadic.h"

/*
    Carries out the finite series evaluation for the logarithm 
    \begin{equation*}
    \sum_{i=1}^{n} a_i y^i
    = \sum_{j=0}^{\ceil{n/b}-1} \Bigl(\sum_{i=1}^b a_{i+jb} y^i\Bigr) y^{jb}
    \end{equation*}
    where $a_i = 1/i$ with the choice $b = \floor{\sqrt{n}}$, 
    all modulo $p^N$, where also $P = p^N$.

    Assumes that $y$ is reduced modulo $p^N$.

    Assumes that $z$ has space for $2d - 1$ coefficients, but 
    sets only the first $d$ to meaningful values on exit.

    Supports aliasing between $y$ and $z$.
 */
static void 
_qadic_log_rectangular_series(fmpz *z, const fmpz *y, slong len, slong n, 
                       const fmpz *a, const slong *j, slong lena, 
                       const fmpz_t p, slong N, const fmpz_t pN)
{
    const slong d = j[lena - 1];

    if (n <= 2)
    {
        if (n == 1)  /* n == 1;  z = y */
        {
            _fmpz_vec_set(z, y, len);
            _fmpz_vec_zero(z + len, d - len);
        }
        else  /* n == 2;  z = y + y^2/2 */
        {
            slong i;
            fmpz *t;

            t = _fmpz_vec_init(2 * len - 1);

            _fmpz_poly_sqr(t, y, len);
            for (i = 0; i < 2 * len - 1; i++)
                if (fmpz_is_even(t + i))
                {
                    fmpz_fdiv_q_2exp(t + i, t + i, 1);
                }
                else  /* => p and t(i) are odd */
                {
                    fmpz_add(t + i, t + i, pN);
                    fmpz_fdiv_q_2exp(t + i, t + i, 1);
                }
            _fmpz_mod_poly_reduce(t, 2 * len - 1, a, j, lena, pN);
            _fmpz_mod_poly_add(z, y, len, t, FLINT_MIN(d, 2 * len - 1), pN);

            _fmpz_vec_clear(t, 2 * len - 1);
        }
    }
    else  /* n >= 3 */
    {
        const slong b = n_sqrt(n);
        const slong k = fmpz_fits_si(p) ? n_flog(n, fmpz_get_si(p)) : 0;

        slong i, h;
        fmpz_t f, pNk;
        fmpz *c, *t, *ypow;

        c    = _fmpz_vec_init(d);
        t    = _fmpz_vec_init(2 * d - 1);
        ypow = _fmpz_vec_init((b + 1) * d + d - 1);
        fmpz_init(f);
        fmpz_init(pNk);

        fmpz_pow_ui(pNk, p, N + k);

        fmpz_one(ypow);
        _fmpz_vec_set(ypow + d, y, len);
        for (i = 2; i <= b; i++)
        {
            _fmpz_mod_poly_mul(ypow + i * d, ypow + (i - 1) * d, d, y, len, pNk);
            _fmpz_mod_poly_reduce(ypow + i * d, d + len - 1, a, j, lena, pNk);
        }

        _fmpz_vec_zero(z, d);

        for (h = (n + (b - 1)) / b - 1; h >= 0; h--)
        {
            const slong hi = FLINT_MIN(b, n - h*b);
            slong w;

            /* Compute inner sum in c */
            fmpz_rfac_uiui(f, 1 + h*b, hi);

            _fmpz_vec_zero(c, d);
            for (i = 1; i <= hi; i++)
            {
                fmpz_divexact_ui(t, f, i + h*b);
                _fmpz_vec_scalar_addmul_fmpz(c, ypow + i * d, d, t);
            }

            /* Multiply c by p^k f */
            w = fmpz_remove(f, f, p);
            _padic_inv(f, f, p, N + k);
            if (w > k)
            {
                fmpz_pow_ui(t, p, w - k);
                _fmpz_vec_scalar_divexact_fmpz(c, c, d, t);
            }
            else if (w < k)
            {
                fmpz_pow_ui(t, p, k - w);
                _fmpz_vec_scalar_mul_fmpz(c, c, d, t);
            }
            _fmpz_vec_scalar_mul_fmpz(c, c, d, f);

            /* Set z = z y^b + c */
            _fmpz_mod_poly_mul(t, z, d, ypow + b * d, d, pNk);
            _fmpz_mod_poly_reduce(t, 2 * d - 1, a, j, lena, pNk);
            _fmpz_vec_add(z, c, t, d);
            _fmpz_vec_scalar_mod_fmpz(z, z, d, pNk);
        }

        fmpz_pow_ui(f, p, k);
        _fmpz_vec_scalar_divexact_fmpz(z, z, d, f);

        fmpz_clear(f);
        fmpz_clear(pNk);
        _fmpz_vec_clear(c, d);
        _fmpz_vec_clear(t, 2 * d - 1);
        _fmpz_vec_clear(ypow, (b + 1) * d + d - 1);

    }
}

void _qadic_log_rectangular(fmpz *z, const fmpz *y, slong v, slong len, 
                            const fmpz *a, const slong *j, slong lena, 
                            const fmpz_t p, slong N, const fmpz_t pN)
{
    const slong d = j[lena - 1];
    const slong n = _padic_log_bound(v, N, p) - 1;

    _qadic_log_rectangular_series(z, y, len, n, a, j, lena, p, N, pN);
    _fmpz_mod_poly_neg(z, z, d, pN);
}

int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const fmpz *p  = (&ctx->pctx)->p;
    const slong d   = qadic_ctx_degree(ctx);
    const slong N   = qadic_prec(rop);
    const slong len = op->length;

    if (op->val < 0)
    {
        return 0;
    }
    else
    {
        fmpz *x;
        fmpz_t pN;
        int alloc, ans;

        x = _fmpz_vec_init(len + 1);
        alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

        /* Set x := (1 - op) mod p^N */
        fmpz_pow_ui(x + len, p, op->val);
        _fmpz_vec_scalar_mul_fmpz(x, op->coeffs, len, x + len);
        fmpz_sub_ui(x, x, 1);
        _fmpz_vec_neg(x, x, len);
        _fmpz_vec_scalar_mod_fmpz(x, x, len, pN);

        if (_fmpz_vec_is_zero(x, len))
        {
            padic_poly_zero(rop);
            ans = 1;
        }
        else
        {
            const slong v = _fmpz_vec_ord_p(x, len, p);

            if (v >= 2 || (*p != WORD(2) && v >= 1))
            {
                if (v >= N)
                {
                    padic_poly_zero(rop);
                }
                else
                {
                    padic_poly_fit_length(rop, d);

                    _qadic_log_rectangular(rop->coeffs, x, v, len, 
                                           ctx->a, ctx->j, ctx->len, p, N, pN);
                    rop->val = 0;

                    _padic_poly_set_length(rop, d);
                    _padic_poly_normalise(rop);
                    padic_poly_canonicalise(rop, p);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        _fmpz_vec_clear(x, len + 1);
        if (alloc)
            fmpz_clear(pN);
        return ans;
    }
}

