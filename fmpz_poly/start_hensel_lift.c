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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/*
    This function takes the local factors (in \code{local_fac}) and 
    Hensel lifts them until they are known mod $p^N$, where $N > 0$. 
    These lifted factors will be stored (in the same ordering) in 
    \code{lifted_fac}. It is assumed that \code{link}, \code{v}, and 
    \code{w} are initialized arrays \code{fmpz_poly_t}'s with at least 
    $2*r - 2$ entries and that $r \geq 2$.  These are done outside of 
    this function so that you can keep them for restarting Hensel lifting 
    later. The product of local factors must be squarefree.

    The return value is an exponent which must be passed to the function 
    \code{_fmpz_poly_continue_hensel_lift()} as \code{prev_exp} if the 
    Hensel lifting is to be resumed.
 */

long _fmpz_poly_start_hensel_lift(fmpz_poly_factor_t lifted_fac, long *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    const nmod_poly_factor_t local_fac, long N)
{
    const long r = local_fac->num;

    long i, preve;
    fmpz_t p, P, big_P;
    fmpz_poly_t monic_f;

    fmpz_init(p);
    fmpz_init(P);
    fmpz_poly_init(monic_f);

    /* Set P := p, monic_f := monic(f) */
    fmpz_set_ui(p, (local_fac->p[0])->mod.n);
    fmpz_set(P, p);

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
        fmpz_t t, Q;

        fmpz_init(Q);
        fmpz_init(t);
        fmpz_pow_ui(Q, P, N);
        fmpz_mod(t, fmpz_poly_lead(f), Q);

        if (fmpz_invmod(t, t, Q))
        {
            printf("Exception in fmpz_poly_start_hensel_lift.\n");
            abort();
        }

        fmpz_poly_scalar_mul_fmpz(monic_f, f, t);
        fmpz_poly_scalar_mod_fmpz(monic_f, monic_f, Q);
        fmpz_clear(Q);
        fmpz_clear(t);
    }

    fmpz_poly_build_hensel_tree(link, v, w, local_fac);

    {
        long *e, n = FLINT_CLOG2(N) + 1;

        e = malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        for (i--; i > 0; i--)
        {
            fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, 
                p, e[i+1], e[i], 1);
        }
        if (N > 1)
        {
            fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, 
                p, e[i+1], e[i], 0);
        }

        preve = e[i+1];

        free(e);
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
            fmpz_poly_scalar_mod_fmpz(lifted_fac->p + (- link[i] - 1), v[i], P);
            lifted_fac->exp[- link[i] - 1] = 1L; 
        }
    }
    lifted_fac->num = r;

    fmpz_clear(p);
    fmpz_clear(P);
    fmpz_poly_clear(monic_f);

    return preve;
}

