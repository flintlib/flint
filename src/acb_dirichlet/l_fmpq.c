/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb_mat.h"
#include "arb_hypgeom.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

int
acb_dirichlet_l_fmpq_use_afe(ulong q, const fmpq_t s, slong prec)
{
    double cutoff;

    if ((slong) fmpz_bits(fmpq_numref(s)) - (slong) fmpz_bits(fmpq_denref(s)) > 20)
        return 0;

    if (fabs(fmpq_get_d(s)) > 10 + prec * 0.01)
        return 0;

    if (q == 1)
    {
        if (fmpz_is_one(fmpq_denref(s)))
            return 0;

        if (fmpz_is_one(fmpq_numref(s)) && fmpz_equal_si(fmpq_denref(s), 2))
            return prec > 32000;

        return prec > 70000;
    }

    if (fmpq_is_zero(s))
        return 0;

    if (fmpz_cmp_ui(fmpq_denref(s), 2) <= 0)
    {
        if (prec > 30000)
            return 1;

        if (fmpq_is_one(s))
            return q > 1000;

        return q > 50;
    }

    cutoff = 15000.0 / q;

    return prec > cutoff;
}

void
acb_dirichlet_l_fmpq(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong q;

    q = (G == NULL) ? 1 : G->q;

    if (acb_dirichlet_l_fmpq_use_afe(q, s, prec))
    {
        acb_dirichlet_l_fmpq_afe(res, s, G, chi, prec);
        if (acb_is_finite(res))
            return;
    }

    acb_set_fmpq(res, s, prec + 10);
    acb_dirichlet_l(res, res, G, chi, prec);
}
