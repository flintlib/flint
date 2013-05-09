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

    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "fmpz_poly_factor.h"

len_t _fmpz_poly_hensel_continue_lift(fmpz_poly_factor_t lifted_fac, 
    len_t *link, fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    len_t prev, len_t curr, len_t N, const fmpz_t p)
{
    const len_t r = lifted_fac->num;

    len_t i, new_prev;

    fmpz_t P;
    fmpz_poly_t monic_f;

    fmpz_init(P);
    fmpz_pow_ui(P, p, N);
    fmpz_poly_init(monic_f);

    if (fmpz_is_one(fmpz_poly_lead(f)))
    {
        fmpz_poly_set(monic_f, f);
    }
    else if (fmpz_cmp_si(fmpz_poly_lead(f), -1) == 0)
    {
        fmpz_poly_neg(monic_f, f);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_mod(t, fmpz_poly_lead(f), P);

        if (fmpz_invmod(t, t, P) == 0)
        {
            printf("Exception (fmpz_poly_continue_hensel_lift).\n");
            abort();
        }

        fmpz_poly_scalar_mul_fmpz(monic_f, f, t);
        fmpz_poly_scalar_mod_fmpz(monic_f, monic_f, P);
        fmpz_clear(t);
    }

    {
        len_t *e, n = 3 + FLINT_FLOG2(N - prev);

        e = flint_malloc(n * sizeof(len_t));

        for (e[i = 0] = N; e[i] > curr; i++)
            e[i + 1] = (e[i] + 1) / 2;
        e[i]   = curr;
        e[i+1] = prev;

        if (prev < curr)
            fmpz_poly_hensel_lift_tree(link, v, w, monic_f, r, p, e[i+1], e[i], -1);

        for (i--; i > 0; i--)
            fmpz_poly_hensel_lift_tree(link, v, w, monic_f, r, p, e[i+1], e[i], 1);

        fmpz_poly_hensel_lift_tree(link, v, w, monic_f, r, p, e[i+1], e[i], 0);

        new_prev = e[i+1];

        flint_free(e);
    }

    /*
        Now everything is lifted to p^N, we just need to 
        insert the factors into their correct places in lifted_fac.
     */
    fmpz_poly_factor_fit_length(lifted_fac, r);

    for(i = 0; i < 2*r - 2; i++)
    { 
        if (link[i] < 0)
        {
            fmpz_poly_scalar_smod_fmpz(lifted_fac->p + (- link[i] - 1), v[i], P);
            lifted_fac->exp[- link[i] - 1] = 1L; 
        }
    }
    lifted_fac->num = r;

    fmpz_clear(P);
    fmpz_poly_clear(monic_f);

    return new_prev;
}

