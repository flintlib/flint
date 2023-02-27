/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "fmpz.h"

#define FLINT_NUM_TINY_DIVISORS FLINT_BITS

const int FLINT_TINY_DIVISORS_SIZE[FLINT_NUM_TINY_DIVISORS] = {
    0,1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6,4,4,2,8,3,4,4,6,2,8,2,
#if FLINT64
    6,4,4,4,9,2,4,4,8,2,8,2,6,6,4,2,10,3,6,4,6,2,8,4,8,4,4,2,12,2,4,6
#endif 
};

const ulong FLINT_TINY_DIVISORS_LOOKUP[FLINT_NUM_TINY_DIVISORS] = {
    UWORD(0x0),UWORD(0x2),UWORD(0x6),0xaUL,UWORD(0x16),UWORD(0x22),0x4eUL,UWORD(0x82),UWORD(0x116),0x20aUL,
    UWORD(0x426),UWORD(0x802),0x105eUL,UWORD(0x2002),UWORD(0x4086),0x802aUL,UWORD(0x10116),UWORD(0x20002),
    0x4024eUL,UWORD(0x80002),UWORD(0x100436),0x20008aUL,UWORD(0x400806),UWORD(0x800002),
    0x100115eUL,UWORD(0x2000022),UWORD(0x4002006),0x800020aUL,UWORD(0x10004096),UWORD(0x20000002),
    0x4000846eUL,UWORD(0x80000002),
#if FLINT64
    UWORD(0x100010116),0x20000080aUL,UWORD(0x400020006),UWORD(0x8000000a2),0x100004125eUL,
    UWORD(0x2000000002),UWORD(0x4000080006),0x800000200aUL,UWORD(0x10000100536),
    UWORD(0x20000000002),0x400002040ceUL,UWORD(0x80000000002),UWORD(0x100000400816),
    0x20000000822aUL,UWORD(0x400000800006),UWORD(0x800000000002),0x100000101115eUL,
    UWORD(0x2000000000082),UWORD(0x4000002000426),0x800000002000aUL,UWORD(0x10000004002016),
    UWORD(0x20000000000002),0x4000000804024eUL,UWORD(0x80000000000822),
    UWORD(0x100000010004196),0x20000000008000aUL,UWORD(0x400000020000006),
    UWORD(0x800000000000002),0x100000004010947eUL,UWORD(0x2000000000000002),
    UWORD(0x4000000080000006),0x800000000020028aUL
#endif
};


void
_arith_divisors(fmpz *res, slong size, fmpz_factor_t factors)
{
    slong i;
    slong *exp = flint_malloc(sizeof(slong) * factors->num);
    slong *exp_max = flint_malloc(sizeof(slong) * factors->num);
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
_arith_divisors_tiny(fmpz_poly_t res, slong n)
{
    slong size;
    slong i, k;

    size = FLINT_TINY_DIVISORS_SIZE[n];

    fmpz_poly_fit_length(res, size);
    i = 0;
    for (k = 1; k <= n; k++)
    {
        if (FLINT_TINY_DIVISORS_LOOKUP[n] & (UWORD(1) << k))
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
    slong i, size, m;
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
