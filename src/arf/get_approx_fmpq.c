/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "arf.h"

static int
cont_frac_step(fmpz_t r, arf_t next, const arf_t current, slong prec, slong tol_exp)
{
    int res = 0;
    arf_get_fmpz(r, current, ARF_RND_FLOOR);
    arf_sub_fmpz(next, current, r, prec, ARF_RND_NEAR);
    if (arf_cmp_2exp_si(next, tol_exp) < 0)
    {
        res = 1;
    }
    else
    {
        arf_ui_div(next, 1, next, prec, ARF_RND_NEAR);
    }
    return res;
}

static void
cont_frac_get_fmpq(fmpq_t c, fmpz* r_vec, slong nb_steps)
{
    slong k;
    fmpq_zero(c);
    fmpq_add_fmpz(c, c, &r_vec[nb_steps-1]);
    for (k = nb_steps-2; k >= 0; k--)
    {
        fmpq_inv(c, c);
        fmpq_add_fmpz(c, c, &r_vec[k]);
    }
}

int arf_get_approx_fmpq(fmpq_t y, const arf_t x, slong prec)
{
    arf_t z;
    int res = 1;
    int stop = 0;
    slong max_steps = prec / 2;
    fmpz* r_vec;
    slong k;

    arf_init(z);
    r_vec = _fmpz_vec_init(max_steps);

    arf_set(z, x);
    k = 0;
    for (k = 0; (k < max_steps) && !stop; k++)
    {
        stop = cont_frac_step(&r_vec[k], z, z, prec, -prec / 8);
    }

    if (k == max_steps)
    {
        res = 0;
    }
    else
    {
        res = 1;
        cont_frac_get_fmpq(y, r_vec, k);
    }

    arf_clear(z);
    _fmpz_vec_clear(r_vec, max_steps);
    return res;
}
