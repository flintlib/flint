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
    Assumes that $1 \leq v$ or $2 \leq v$ as $p$ is even 
    or odd, respectively, and that $v < N < 2^{f-2}$ where 
    $f$ is \code{FLINT_BITS}.

    Under the assumption that $1 \leq v < N$, or $2 \leq v < N$, 
    one can easily prove that with $c = N - \floor{\log_p v}$ 
    the number $b = \ceil{(c + \ceil{\log_p c} + 1) / b}$ is such 
    that for all $i \geq b$, $i v - \ord_p(i) \geq N$.

    Under the additional condition that $N < 2^{f-2}$ one can 
    show that the code branch for primes that fit into a 
    \code{signed long} does not cause overflow.  Moreover, 
    independently of this, it follows that the above value $b$ 
    is less than $2^{f-1}$.

    In the first branch, we have that $b v - log_p(b) \geq N$.
    We need to show that we can replace $\log_p$ by $\ord_p$ here.  
    That is, we need that $iv - \ord_p(i) \geq iv - \log_p(i) \geq N$, 
    i.e., $\log_p(i) \geq \ord_p(i)$, which is true.  We then work 
    backwards to find the first $i$ such that this fails, then 
    using that the function is strictly increasing for $i \geq 2$.

    In the second branch we use that using signed indices in the 
    summation is still sufficient and hence that all terms $1/i$ 
    are units.
    Then $ord_p(x^i/i) \geq N$ provided that $i v \geq N$.
 */

len_t _padic_log_bound(len_t v, len_t N, const fmpz_t prime)
{
    if (N >= (1L << (FLINT_BITS - 2)))
    {
        printf("Exception (_padic_log_bound).  N = %ld is too large.\n", N);
        abort();
    }

    if (fmpz_fits_si(prime))
    {
        len_t b, c, p = fmpz_get_si(prime);

        c = N - n_flog(v, p);
        b = ((c + n_clog(c, p) + 1) + (v - 1)) / v;

        while (--b >= 2)
        {
            len_t t = b * v - n_clog(b, p);

            if (t < N)
                return b + 1;
        }

        return 2;
    }
    else
    {
        return (N + v - 1) / v;
    }
}

void _padic_log(fmpz_t z, const fmpz_t y, len_t v, const fmpz_t p, len_t N)
{
    if (N < (1L << 9) / (len_t) fmpz_bits(p))
    {
        _padic_log_rectangular(z, y, v, p, N);
    }
    else
    {
        _padic_log_balanced(z, y, v, p, N);
    }
}

int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
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
                    _padic_log(padic_unit(rop), x, v, p, N);
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

