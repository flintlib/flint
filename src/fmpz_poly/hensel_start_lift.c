/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

slong _fmpz_poly_hensel_start_lift(fmpz_poly_factor_t lifted_fac, slong *link,
    fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f,
    const nmod_poly_factor_t local_fac, slong N)
{
    const slong r = local_fac->num;

    slong i, preve;
    fmpz_t p, P;
    fmpz_poly_t monic_f;

    fmpz_init(p);
    fmpz_init(P);
    fmpz_poly_init(monic_f);

    fmpz_set_ui(p, (local_fac->p + 0)->mod.n);
    fmpz_pow_ui(P, p, N);

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
            flint_throw(FLINT_ERROR, "Exception (fmpz_poly_start_hensel_lift).\n");
        }

        fmpz_poly_scalar_mul_fmpz(monic_f, f, t);
        fmpz_poly_scalar_mod_fmpz(monic_f, monic_f, P);
        fmpz_clear(t);
    }

    fmpz_poly_hensel_build_tree(link, v, w, local_fac);

    {
        slong *e, n = FLINT_CLOG2(N) + 1;

        e = flint_malloc(n * sizeof(slong));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        for (i--; i > 0; i--)
        {
            fmpz_poly_hensel_lift_tree(link, v, w, monic_f, r,
                p, e[i+1], e[i], 1);
        }
        if (N > 1)
        {
            fmpz_poly_hensel_lift_tree(link, v, w, monic_f, r,
                p, e[i+1], e[i], 0);
        }

        preve = e[i+1];

        flint_free(e);
    }

    /*
        Now everything is lifted to p^N, we just need to
        insert the factors into their correct places in lifted_fac.
     */
    fmpz_poly_factor_fit_length(lifted_fac, r);

    for (i = 0; i < 2*r - 2; i++)
    {
        if (link[i] < 0)
        {
            fmpz_poly_scalar_smod_fmpz(lifted_fac->p + (- link[i] - 1), v[i], P);
            lifted_fac->exp[- link[i] - 1] = WORD(1);
        }
    }
    lifted_fac->num = r;

    fmpz_clear(p);
    fmpz_clear(P);
    fmpz_poly_clear(monic_f);

    return preve;
}

