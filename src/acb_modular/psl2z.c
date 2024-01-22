/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

void
psl2z_inv(psl2z_t h, const psl2z_t g)
{
    if (h != g)
        psl2z_set(h, g);

    if (fmpz_is_zero(&h->c) && fmpz_sgn(&h->a) > 0)
    {
        fmpz_inplace_neg(&h->b);
        fmpz_swap(&h->d, &h->a);
    }
    else
    {
        fmpz_swap(&h->a, &h->d);
        fmpz_inplace_neg(&h->a);
        fmpz_inplace_neg(&h->d);
    }
}

int
psl2z_is_correct(const psl2z_t g)
{
    int res;
    fmpz_t t;

    if (fmpz_sgn(&g->c) < 0)
        return 0;

    if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) <= 0)
        return 0;

    fmpz_init(t);
    fmpz_mul(t, &g->a, &g->d);
    fmpz_submul(t, &g->b, &g->c);
    res = fmpz_is_one(t);
    fmpz_clear(t);
    return res;
}

int
psl2z_is_one(const psl2z_t g)
{
    return fmpz_is_one(&g->a) && fmpz_is_zero(&g->b) &&
            fmpz_is_zero(&g->c) && fmpz_is_one(&g->d);
}

void
psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g)
{
    if (h == f || h == g)
    {
        psl2z_t t;
        psl2z_init(t);
        psl2z_mul(t, f, g);
        psl2z_swap(t, h);
        psl2z_clear(t);
        return;
    }

    fmpz_mul(&h->a, &f->a, &g->a);
    fmpz_addmul(&h->a, &f->b, &g->c);

    fmpz_mul(&h->b, &f->a, &g->b);
    fmpz_addmul(&h->b, &f->b, &g->d);

    fmpz_mul(&h->c, &f->c, &g->a);
    fmpz_addmul(&h->c, &f->d, &g->c);

    fmpz_mul(&h->d, &f->c, &g->b);
    fmpz_addmul(&h->d, &f->d, &g->d);

    if (fmpz_sgn(&h->c) < 0 || (fmpz_is_zero(&h->c) && fmpz_sgn(&h->d) < 0))
    {
        fmpz_inplace_neg(&h->a);
        fmpz_inplace_neg(&h->b);
        fmpz_inplace_neg(&h->c);
        fmpz_inplace_neg(&h->d);
    }
}

void
psl2z_randtest(psl2z_t g, flint_rand_t state, slong bits)
{
    bits = FLINT_MAX(bits, 1);

    fmpz_randtest(&g->a, state, bits);
    fmpz_randtest(&g->b, state, bits);

    if (fmpz_is_zero(&g->a) && fmpz_is_zero(&g->b))
    {
        psl2z_one(g);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_xgcd(t, &g->d, &g->c, &g->a, &g->b);
        fmpz_divexact(&g->a, &g->a, t);
        fmpz_divexact(&g->b, &g->b, t);

        if (fmpz_sgn(&g->c) < 0)
            fmpz_inplace_neg(&g->c);
        else
            fmpz_inplace_neg(&g->b);

        if (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) < 0)
        {
            fmpz_inplace_neg(&g->a);
            fmpz_inplace_neg(&g->b);
            fmpz_inplace_neg(&g->c);
            fmpz_inplace_neg(&g->d);
        }

        fmpz_clear(t);
    }
}
