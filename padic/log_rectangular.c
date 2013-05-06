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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"

/*
    Carries out the finite series evaluation for the logarithm 
    \begin{equation*}
    \sum_{i=1}^{n} a_i x^i
    = \sum_{j=0}^{\ceil{n/b} - 1} \Bigl( \sum_{i=1}^b a_{i+jb} x^i \Bigr) x^{jb}
    \end{equation*}
    where $a_i = 1/i$ with the choice $b = \floor{\sqrt{n}}$, 
    all modulo $p^N$, where also $P = p^N$.

    Does not support aliasing.
 */
static void 
_padic_log_rectangular_series(fmpz_t z, const fmpz_t y, long n, 
                              const fmpz_t p, long N, const fmpz_t P0)
{
    if (n <= 2)
    {
        if (n == 1)
        {
            fmpz_mod(z, y, P0);
        }
        else  /* n == 2;  z = y(1 + y/2) */
        {
            if (fmpz_is_even(y))
            {
                fmpz_fdiv_q_2exp(z, y, 1);
            }
            else  /* => p and y are odd */
            {
                fmpz_add(z, y, P0);
                fmpz_fdiv_q_2exp(z, z, 1);
            }
            fmpz_add_ui(z, z, 1);
            fmpz_mul(z, z, y);
            fmpz_mod(z, z, P0);
        }
    }
    else
    {
        const long b = n_sqrt(n);
        const long k = fmpz_fits_si(p) ? n_flog(n, fmpz_get_si(p)) : 0;

        long i, j;
        fmpz_t c, f, t, P1;
        fmpz *ypow;

        ypow = _fmpz_vec_init(b + 1);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(t);
        fmpz_init(P1);

        fmpz_pow_ui(P1, p, N + k);

        fmpz_one(ypow + 0);
        for (i = 1; i <= b; i++)
        {
            fmpz_mul(ypow + i, ypow + (i - 1), y);
            fmpz_mod(ypow + i, ypow + i, P1);
        }

        fmpz_zero(z);

        for (j = (n + (b - 1)) / b - 1; j >= 0; j--)
        {
            const long hi = FLINT_MIN(b, n - j*b);
            long w;

            /* Compute inner sum in c */
            fmpz_rfac_uiui(f, 1 + j*b, hi);

            fmpz_zero(c);
            for (i = 1; i <= hi; i++)
            {
                fmpz_divexact_ui(t, f, i + j*b);
                fmpz_addmul(c, t, ypow + i);
            }

            /* Multiply c by p^k f^{-1} */
            w = fmpz_remove(f, f, p);
            _padic_inv(f, f, p, N);
            if (w > k)
            {
                fmpz_pow_ui(t, p, w - k);
                fmpz_divexact(c, c, t);
            }
            else
            {
                fmpz_pow_ui(t, p, k - w);
                fmpz_mul(c, c, t);
            }
            fmpz_mul(c, c, f);

            /* Set z = z y^b + c */
            fmpz_mul(t, z, ypow + b);
            fmpz_add(z, c, t);
            fmpz_mod(z, z, P1);
        }

        fmpz_pow_ui(f, p, k);
        fmpz_divexact(z, z, f);

        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(t);
        fmpz_clear(P1);
        _fmpz_vec_clear(ypow, b + 1);
    }
}

void _padic_log_rectangular(fmpz_t z, const fmpz_t y, long v, const fmpz_t p, long N)
{
    const long n = _padic_log_bound(v, N, p) - 1;
    fmpz_t pN;

    fmpz_init(pN);
    fmpz_pow_ui(pN, p, N);

    _padic_log_rectangular_series(z, y, n, p, N, pN);

    fmpz_sub(z, pN, z);
    fmpz_clear(pN);
}

int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const fmpz *p = ctx->p;
    const long N  = padic_prec(rop);

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
                    _padic_log_rectangular(padic_unit(rop), x, v, p, N);
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

