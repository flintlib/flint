/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

/*
    TODO:  Move this bit of code into "padic".
 */
static void __padic_reduce(fmpz_t u, slong *v, slong N, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(u))
    {
        if (*v < N)
        {
            int alloc;
            fmpz_t pow;

            alloc = _padic_ctx_pow_ui(pow, N - *v, ctx);
            fmpz_mod(u, u, pow);
            if (alloc)
                fmpz_clear(pow);
        }
        else
        {
            fmpz_zero(u);
            *v = 0;
        }
    }
}

/*
    Evaluates the polynomial $F(x) = p^w f(x)$ at $x = p^b a$, setting 
    $y = p^v u$ to the result reduced modulo $p^N$.

    Suppose first that $b \geq 0$, in which case we can quickly relay 
    the call to the \code{fmpz_mod_poly} module.  Namely, we need to 
    compute $f(x) \bmod {p^{N-w}}$ where we know $f(x)$ to be integral.

    Otherwise, suppose now that $b < 0$ and we still wish to evaluate 
    $f(x) \bmod {p^{N-w}}$.
    \begin{align*}
    f(x) & = \sum_{i = 0}^{n} a_i x^i \\
         & = \sum_{i = 0}^{n} a_i p^{i b} a^i \\
    \intertext{Multiplying through by $p^{- n b} \in \mathbf{Z}$, }
    p^{-nb} f(x) & = \sum_{i = 0}^{n} a_i p^{-(n-i)b} a
    \end{align*}
    which leaves the right hand side integral.  As we want to 
    compute $f(x)$ to precision $N-w$, we have to compute 
    $p^{-nb} f(x)$ to precision $N-w-nb$.
 */

void _padic_poly_evaluate_padic(fmpz_t u, slong *v, slong N,  
                                const fmpz *poly, slong val, slong len, 
                                const fmpz_t a, slong b, const padic_ctx_t ctx)
{
    if (len == 0)
    {
        fmpz_zero(u);
        *v = 0;
    }
    else if (len == 1)
    {
        fmpz_set(u, poly);
        *v = val;

        __padic_reduce(u, v, N, ctx);
    }
    else if (b >= 0)
    {
        if (val >= N)
        {
            fmpz_zero(u);
            *v = 0;
        }
        else 
        {
            fmpz_t x;
            fmpz_t pow;
            int alloc;

            fmpz_init(x);
            alloc = _padic_ctx_pow_ui(pow, N - val, ctx);

            fmpz_pow_ui(x, ctx->p, b);
            fmpz_mul(x, x, a);

            _fmpz_mod_poly_evaluate_fmpz(u, poly, len, x, pow);
            if (!fmpz_is_zero(u))
                *v = val + _fmpz_remove(u, ctx->p, ctx->pinv);
            else
                *v = 0;

            fmpz_clear(x);
            if (alloc)
                fmpz_clear(pow);
        }
    }
    else  /* b < 0 */
    {
        const slong n = len - 1;

        if (val + n*b >= N)
        {
            fmpz_zero(u);
            *v = 0;
        }
        else
        {
            fmpz_t pow;
            int alloc;
            slong i;
            fmpz_t s, t;

            fmpz *vec = _fmpz_vec_init(len);
            fmpz_init(s);
            fmpz_init(t);

            alloc = _padic_ctx_pow_ui(pow, N - val - n*b, ctx);

            fmpz_pow_ui(s, ctx->p, -b);
            fmpz_one(t);
            fmpz_set(vec + (len - 1), poly + (len - 1));
            for (i = len - 2; i >= 0; i--)
            {
                fmpz_mul(t, t, s);
                fmpz_mul(vec + i, poly + i, t);
            }

            _fmpz_mod_poly_evaluate_fmpz(u, vec, len, a, pow);
            if (!fmpz_is_zero(u))
                *v = val + n*b + _fmpz_remove(u, ctx->p, ctx->pinv);
            else
                *v = 0;

            if (alloc)
                fmpz_clear(pow);
            fmpz_clear(s);
            fmpz_clear(t);
            _fmpz_vec_clear(vec, len);
        }
    }
}

void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, 
                                          const padic_t x, const padic_ctx_t ctx)
{
    if (y == x)
    {
        padic_t t;

        padic_init2(t, padic_prec(y));
        _padic_poly_evaluate_padic(padic_unit(t), &padic_val(t), padic_prec(t), 
                                   poly->coeffs, poly->val, poly->length, 
                                   padic_unit(x), padic_val(x), ctx);
        padic_swap(y, t);
        padic_clear(t);
    }
    else
    {
        _padic_poly_evaluate_padic(padic_unit(y), &padic_val(y), padic_prec(y), 
                                   poly->coeffs, poly->val, poly->length, 
                                   padic_unit(x), padic_val(x), ctx);
    }
}

