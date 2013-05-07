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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "arith.h"
#include "ulong_extras.h"

#define FLINT_NUM_TINY_DIVISORS FLINT_BITS

const int FLINT_TINY_DIVISORS_SIZE[FLINT_NUM_TINY_DIVISORS] = {
    0,1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6,4,4,2,8,3,4,4,6,2,8,2,
#if FLINT64
    6,4,4,4,9,2,4,4,8,2,8,2,6,6,4,2,10,3,6,4,6,2,8,4,8,4,4,2,12,2,4,6
#endif 
};

const ulong FLINT_TINY_DIVISORS_LOOKUP[FLINT_NUM_TINY_DIVISORS] = {
    0x0UL,0x2UL,0x6UL,0xaUL,0x16UL,0x22UL,0x4eUL,0x82UL,0x116UL,0x20aUL,
    0x426UL,0x802UL,0x105eUL,0x2002UL,0x4086UL,0x802aUL,0x10116UL,0x20002UL,
    0x4024eUL,0x80002UL,0x100436UL,0x20008aUL,0x400806UL,0x800002UL,
    0x100115eUL,0x2000022UL,0x4002006UL,0x800020aUL,0x10004096UL,0x20000002UL,
    0x4000846eUL,0x80000002UL,
#if FLINT64
    0x100010116UL,0x20000080aUL,0x400020006UL,0x8000000a2UL,0x100004125eUL,
    0x2000000002UL,0x4000080006UL,0x800000200aUL,0x10000100536UL,
    0x20000000002UL,0x400002040ceUL,0x80000000002UL,0x100000400816UL,
    0x20000000822aUL,0x400000800006UL,0x800000000002UL,0x100000101115eUL,
    0x2000000000082UL,0x4000002000426UL,0x800000002000aUL,0x10000004002016UL,
    0x20000000000002UL,0x4000000804024eUL,0x80000000000822UL,
    0x100000010004196UL,0x20000000008000aUL,0x400000020000006UL,
    0x800000000000002UL,0x100000004010947eUL,0x2000000000000002UL,
    0x4000000080000006UL,0x800000000020028aUL
#endif
};


void
_arith_divisors(fmpz *res, long size, fmpz_factor_t factors)
{
    long i;
    long *exp = flint_malloc(sizeof(long) * factors->num);
    long *exp_max = flint_malloc(sizeof(long) * factors->num);
    fmpz *powers = _fmpz_vec_init(factors->num);
    fmpz_t d;

    for (i = 0; i < factors->num; i++)
    {
        exp[i] = 0;
        fmpz_set(powers + i, factors->p + i);
        exp_max[i] = factors->exp[i];
        fmpz_pow_ui(powers + i, powers + i, exp_max[i]);
    }

    fmpz_init(d);
    fmpz_one(res);
    fmpz_one(d);
    res++;

    i = 0;
    while (1)
    {
        while (1)
        {
            if (i == factors->num)
                goto all_done;
            if (exp[i] < exp_max[i])
            {
                exp[i]++;
                fmpz_mul(d, d, factors->p + i);
                i = 0;
                break;
            }
            else
            {
                exp[i] = 0;
                fmpz_divexact(d, d, powers+i);
                i += 1;
            }
        }
        fmpz_set(res, d);
        res++;
    }

    all_done:
    fmpz_clear(d);
    flint_free(exp);
    flint_free(exp_max);
    _fmpz_vec_clear(powers, factors->num);
}


void
_arith_divisors_tiny(fmpz_poly_t res, long n)
{
    long size;
    long i, k;

    size = FLINT_TINY_DIVISORS_SIZE[n];

    fmpz_poly_fit_length(res, size);
    i = 0;
    for (k = 1; k <= n; k++)
    {
        if (FLINT_TINY_DIVISORS_LOOKUP[n] & (1UL << k))
        {
            fmpz_poly_set_coeff_si(res, i, k);
            i++;
        }
    }
    _fmpz_poly_set_length(res, size);
    return;
}

void
arith_divisors(fmpz_poly_t res, const fmpz_t n)
{
    long i, size, m;
    fmpz_factor_t factors;

    if (!COEFF_IS_MPZ(*n))
    {
        m = fmpz_get_si(n);
        if (-FLINT_NUM_TINY_DIVISORS < m && m < FLINT_NUM_TINY_DIVISORS)
        {
            _arith_divisors_tiny(res, FLINT_ABS(m));
            return;
        }
    }

    fmpz_factor_init(factors);
    fmpz_factor(factors, n);

    /* TODO: check for overflow for huge n */
    size = 1;
    for (i = 0; i < factors->num; i++)
        size *= factors->exp[i] + 1;

    fmpz_poly_fit_length(res, size);
    _arith_divisors(res->coeffs, size, factors);
    _fmpz_poly_set_length(res, size);
    _fmpz_vec_sort(res->coeffs, size);

    fmpz_factor_clear(factors);
}
