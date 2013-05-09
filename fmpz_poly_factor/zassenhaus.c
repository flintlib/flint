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

#include <stdlib.h>
#include "fmpz_poly_factor.h"

#define TRACE_ZASSENHAUS 0

/*
    Let $f$ be a polynomial of degree $m = \deg(f) \geq 2$. 
    If another polynomial $g$ divides $f$ then, for all 
    $0 \leq j \leq \deg(g)$, 
    \begin{equation*}
    \abs{b_j} \leq \binom{n-1}{j} \abs{f} + \binom{n-1}{j-1} \abs{a_m}
    \end{equation*}
    where $\abs{f}$ denotes the $2$-norm of $f$.  This bound 
    is due to Mignotte, see e.g., Cohen p.\ 134.

    This function sets $B$ such that, for all $0 \leq j \leq \deg(g)$, 
    $\abs{b_j} \leq B$.

    Consequently, when proceeding with Hensel lifting, we 
    proceed to choose an $a$ such that $p^a \geq 2 B + 1$, 
    e.g., $a = \ceil{\log_p(2B + 1)}$.

    Note that the formula degenerates for $j = 0$ and $j = n$ 
    and so in this case we use that the leading (resp.\ constant) 
    term of $g$ divides the leading (resp.\ constant) term of $f$.
 */
static void _fmpz_poly_factor_mignotte(fmpz_t B, const fmpz *f, len_t m)
{
    len_t j;
    fmpz_t b, f2, lc, s, t;

    fmpz_init(b);
    fmpz_init(f2);
    fmpz_init(lc);
    fmpz_init(s);
    fmpz_init(t);

    for (j = 0; j <= m; j++)
        fmpz_addmul(f2, f + j, f + j);
    fmpz_sqrt(f2, f2);
    fmpz_add_ui(f2, f2, 1);

    fmpz_abs(lc, f + m);

    fmpz_abs(B, f + 0);

    /*  We have $b = \binom{m-1}{j-1}$ on loop entry and 
        $b = \binom{m-1}{j}$ on exit. */
    fmpz_set_ui(b, m-1);
    for (j = 1; j < m; j++)
    {
        fmpz_mul(t, b, lc);

        fmpz_mul_ui(b, b, m - j);
        fmpz_divexact_ui(b, b, j);

        fmpz_mul(s, b, f2);
        fmpz_add(s, s, t);
        if (fmpz_cmp(B, s) < 0)
            fmpz_set(B, s);
    }

    if (fmpz_cmp(B, lc) < 0)
        fmpz_set(B, lc);

    fmpz_clear(b);
    fmpz_clear(f2);
    fmpz_clear(lc);
    fmpz_clear(s);
    fmpz_clear(t);
}

static void fmpz_poly_factor_mignotte(fmpz_t B, const fmpz_poly_t f)
{
    _fmpz_poly_factor_mignotte(B, f->coeffs, f->length - 1);
}

void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, 
                                  len_t exp, const fmpz_poly_t f, len_t cutoff)
{
    const len_t lenF = f->length;

    #if TRACE_ZASSENHAUS == 1
    printf("\n[Zassenhaus]\n");
    printf("|f = "), fmpz_poly_print(f), printf("\n");
    #endif

    if (lenF == 2)
    {
        fmpz_poly_factor_insert(final_fac, f, exp);
    }
    else
    {
        len_t i;
        len_t r = lenF;
        mp_limb_t p = 2;
        nmod_poly_t d, g, t;
        nmod_poly_factor_t fac;

        nmod_poly_factor_init(fac);
        nmod_poly_init_preinv(t, 0, 0);
        nmod_poly_init_preinv(d, 0, 0);
        nmod_poly_init_preinv(g, 0, 0);

        for (i = 0; i < 3; i++)
        {
            for ( ; ; p = n_nextprime(p, 0))
            {
                nmod_t mod;

                nmod_init(&mod, p);
                d->mod = mod;
                g->mod = mod;
                t->mod = mod;

                fmpz_poly_get_nmod_poly(t, f);
                if (t->length == lenF)
                {
                    nmod_poly_derivative(d, t);
                    nmod_poly_gcd(g, t, d);

                    if (nmod_poly_is_one(g))
                    {
                        nmod_poly_factor_t temp_fac;

                        nmod_poly_factor_init(temp_fac);
                        nmod_poly_factor(temp_fac, t);

                        if (temp_fac->num <= r)
                        {
                            r = temp_fac->num;
                            nmod_poly_factor_set(fac, temp_fac);
                        }
                        nmod_poly_factor_clear(temp_fac);
                        break;
                    }
                }
            }
            p = n_nextprime(p, 0);
        }
        nmod_poly_clear(d);
        nmod_poly_clear(g);
        nmod_poly_clear(t);

        if (r > cutoff)
        {
            printf("Exception (fmpz_poly_factor_zassenhaus). r > cutoff.\n");
            nmod_poly_factor_clear(fac);
            abort();
        }
        else if (r == 1)
        {
            fmpz_poly_factor_insert(final_fac, f, exp);
        }
        else
        {
            len_t a;
            fmpz_poly_factor_t lifted_fac;
            fmpz_poly_factor_init(lifted_fac);

            p = (fac->p + 0)->mod.n;
            {
                fmpz_t B;
                fmpz_init(B);
                fmpz_poly_factor_mignotte(B, f);
                fmpz_mul_ui(B, B, 2);
                fmpz_add_ui(B, B, 1);
                a = fmpz_clog_ui(B, p);
                fmpz_clear(B);
            }

            /* TODO: Check if use_Hoeij_Novocin and try smaller a. */
            fmpz_poly_hensel_lift_once(lifted_fac, f, fac, a);

            #if TRACE_ZASSENHAUS == 1
            printf("|p = %ld, a = %ld\n", p, a);
            printf("|Pre hensel lift factorisation (nmod_poly):\n");
            nmod_poly_factor_print(fac);
            printf("|Post hensel lift factorisation (fmpz_poly):\n");
            fmpz_poly_factor_print(lifted_fac);
            #endif

            /* Recombination */
            {
                fmpz_t P;
                fmpz_init(P);
                fmpz_set_ui(P, p);
                fmpz_pow_ui(P, P, a);

                fmpz_poly_factor_zassenhaus_recombination(final_fac, lifted_fac, f, P, exp);

                fmpz_clear(P);
            }

            fmpz_poly_factor_clear(lifted_fac);
        }
        nmod_poly_factor_clear(fac);
    }
}

void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    const len_t lenG = G->length;
    fmpz_poly_t g;

    if (lenG == 0)
    {
        fmpz_set_ui(&fac->c, 0);
        return;
    }
    if (lenG == 1)
    {
        fmpz_set(&fac->c, G->coeffs);
        return;
    }

    fmpz_poly_init(g);

    if (lenG == 2)
    {
        fmpz_poly_content(&fac->c, G);
        if (fmpz_sgn(fmpz_poly_lead(G)) < 0)
            fmpz_neg(&fac->c, &fac->c);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);
        fmpz_poly_factor_insert(fac, g, 1);
    }
    else
    {
        len_t j, k;
        fmpz_poly_factor_t sq_fr_fac;

        /* Does a presearch for a factor of form x^k */
        for (k = 0; fmpz_is_zero(G->coeffs + k); k++) ;

        if (k != 0)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_set_coeff_ui(t, 1, 1);
            fmpz_poly_factor_insert(fac, t, k);
            fmpz_poly_clear(t);
        }

        fmpz_poly_shift_right(g, G, k);

        /* Could make other tests for x-1 or simple things 
           maybe take advantage of the composition algorithm */
        fmpz_poly_factor_init(sq_fr_fac);
        fmpz_poly_factor_squarefree(sq_fr_fac, g);

        fmpz_set(&fac->c, &sq_fr_fac->c);

        /* Factor each square-free part */
        for (j = 0; j < sq_fr_fac->num; j++)
            _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j], sq_fr_fac->p + j, 10);

        fmpz_poly_factor_clear(sq_fr_fac);
    }
    fmpz_poly_clear(g);
}

#undef TRACE_ZASSENHAUS

