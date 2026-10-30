/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb.h"
#include "acb_dirichlet.h"

/*
    The approximate functional equation is implemented in l_fmpq_afe_jet.c,
    which handles power series in s; the value is the case len = 1.
*/
void
acb_dirichlet_l_fmpq_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong q;
    int parity;

    q = (G == NULL) ? 1 : G->q;
    parity = (G == NULL) ? 0 : dirichlet_parity_char(G, chi);

    /* Division by gamma((s+parity)/2) at a pole. */
    if (fmpz_is_one(fmpq_denref(s)))
    {
        const fmpz * n = fmpq_numref(s);

        if ((parity == 0 && fmpz_sgn(n) <= 0 && fmpz_is_even(n)) ||
            (parity == 1 && fmpz_sgn(n) < 0 && fmpz_is_odd(n)))
        {
            /* Special case for zeta (with a principal character modulo
               q > 1, the Euler factor at p | q vanishes) */
            if (q == 1 && fmpz_is_zero(n))
                acb_set_d(res, -0.5);
            else
                acb_zero(res);
            return;
        }
    }

    acb_dirichlet_l_fmpq_afe_jet(res, s, G, chi, 0, 1, prec);
}
