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
#include "ulong_extras.h"

/*
    This function restarts a stopped Hensel lift.

    It lifts from \code{curr} to $N$. It also requires \code{prev} 
    (to lift the inverses) given as the return value of the function 
    \code{_fmpz_poly_start_hensel_lift()} or the function 
    \code{_fmpz_poly_continue_hensel_lift()}. The curr lifted factors 
    are supplied in \code{lifted_fac} and upon return are updated
    there. As usual \code{link}, \code{v}, and \code{w} describe the 
    curr Hensel tree, $r$ is the number of local factors and $p$ is 
    the small prime modulo whose power we are lifting to. It is required 
    that \code{curr} be at least $1$ and that \code{N > curr}.
 */

long _fmpz_poly_continue_hensel_lift(fmpz_poly_factor_t lifted_fac, 
    long *link, fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    long prev, long curr, long N, const fmpz_t p)
{
    const long r = lifted_fac->num;

    long i, new_prev;

    fmpz_t P;
    fmpz_poly_t monic_f;

    fmpz_init_set(P, p);
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
        fmpz_t t, Q;

        fmpz_init(Q);
        fmpz_init(t);
        fmpz_pow_ui(Q, P, N);
        fmpz_mod(t, fmpz_poly_lead(f), Q);

        if (fmpz_invmod(t, t, Q))
        {
            printf("Exception in fmpz_poly_continue_hensel_lift.\n");
            abort();
        }

        fmpz_poly_scalar_mul_fmpz(monic_f, f, t);
        fmpz_poly_scalar_mod_fmpz(monic_f, monic_f, Q);
        fmpz_clear(Q);
        fmpz_clear(t);
    }

    {
        long *e, n = 2 + FLINT_FLOG2(N - prev);

        e = malloc(n * sizeof(long));

        for (e[i = 0] = N; e[i] > curr; i++)
            e[i + 1] = (e[i] + 1) / 2;
        e[i]   = curr;
        e[i+1] = prev;

        fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, p, e[i+1], e[i], -1);

        for (i--; i > 0; i--)
            fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, p, e[i+1], e[i], 1);

        fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, p, e[i+1], e[i], 0);

        new_prev = e[i+1];

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

    for(i = 0; i < 2*r - 2; i++)
    { 
        if (link[i] < 0)
        {
            fmpz_poly_scalar_mod_fmpz(lifted_fac->p + (- link[i] - 1), v[i], P);
            lifted_fac->exp[- link[i] - 1] = 1L; 
        }
    }
    lifted_fac->num = r;

    fmpz_clear(P);
    fmpz_poly_clear(monic_f);

    return new_prev;
}

