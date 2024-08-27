/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_ql_recombine(acb_ptr th, acb_srcptr th0, acb_srcptr cofactors,
    const slong * pts, slong nb, const arf_t err, slong fullprec,
    slong s, ulong a, int all, slong g, slong prec)
{
    slong n = 1 << g;
    slong n0 = 1 << s;
    slong nba = 1 << (g - s);
    slong nbth0 = (all ? n0 * n0 : n0);
    acb_ptr aux;
    acb_t x;
    ulong a0, b0, b1, b;
    slong k;

    acb_init(x);
    aux = _acb_vec_init(nbth0);

    /* Sum contributions of each point weighted by cofactor */
    for (k = 0; k < nb; k++)
    {
        _acb_vec_scalar_mul(aux, th0 + k * nbth0, nbth0, &cofactors[k], prec);

        for (a0 = 0; a0 < n0; a0++)
        {
            if (all)
            {
                /* writing ab = a0 a1 b0 b1, we modify all entries with
                   a1 = a, with sign i^(b1 . pt) */
                for (b = 0; b < n; b++)
                {
                    b1 = b % nba;
                    b0 = b >> (g - s);
                    acb_mul_i_pow_si(x, &aux[a0 * n0 + b0],
                        2 * acb_theta_char_dot_slong(b1, pts + k * (g - s), g - s)
                        + acb_theta_char_dot(b1, a, g - s));
                    acb_add(&th[(a0 << (g + g - s)) + (a << g) + b],
                        &th[(a0 << (g + g - s)) + (a << g) + b], x, fullprec);
                }
            }
            else
            {
                acb_add(&th[(a0 << (g - s)) + a], &th[(a0 << (g - s)) + a],
                    &aux[a0], fullprec);
            }
        }
    }
    /* Add error */
    for (a0 = 0; a0 < n0; a0++)
    {
        if (all)
        {
            for (b = 0; b < n; b++)
            {
                acb_add_error_arf(&th[(a0 << (g + g - s)) + (a << g) + b], err);
            }
        }
        else
        {
            acb_add_error_arf(&th[(a0 << (g - s)) + a], err);
        }
    }

    acb_clear(x);
    _acb_vec_clear(aux, nbth0);
}
