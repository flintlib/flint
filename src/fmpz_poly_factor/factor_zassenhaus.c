/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"

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
void _fmpz_poly_factor_mignotte(fmpz_t B, const fmpz *f, slong m)
{
    slong j;
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

void fmpz_poly_factor_mignotte(fmpz_t B, const fmpz_poly_t f)
{
    _fmpz_poly_factor_mignotte(B, f->coeffs, f->length - 1);
}

void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, 
               slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij)
{
    const slong lenF = f->length;

    #if TRACE_ZASSENHAUS == 1
    flint_printf("\n[Zassenhaus]\n");
    flint_printf("|f = "), fmpz_poly_print(f), flint_printf("\n");
    #endif

    if (lenF < 5)
    {
        if (lenF < 3)
            fmpz_poly_factor_insert(final_fac, f, exp);
        else if (lenF == 3)
            _fmpz_poly_factor_quadratic(final_fac, f, exp);
        else
            _fmpz_poly_factor_cubic(final_fac, f, exp);

        return;
    }
    else
    {
        slong i, j;
        slong r = lenF;
        mp_limb_t p = 2;
        nmod_poly_t d, g, t;
        nmod_poly_factor_t fac;
        zassenhaus_prune_t Z;

        zassenhaus_prune_init(Z);
        nmod_poly_factor_init(fac);
        nmod_poly_init_preinv(t, 1, 0);
        nmod_poly_init_preinv(d, 1, 0);
        nmod_poly_init_preinv(g, 1, 0);

        zassenhaus_prune_set_degree(Z, lenF - 1);

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
                if (t->length == lenF && t->coeffs[0] != 0)
                {
                    nmod_poly_derivative(d, t);
                    nmod_poly_gcd(g, t, d);

                    if (nmod_poly_is_one(g))
                    {
                        nmod_poly_factor_t temp_fac;

                        nmod_poly_factor_init(temp_fac);
                        nmod_poly_factor(temp_fac, t);

                        zassenhaus_prune_start_add_factors(Z);
                        for (j = 0; j < temp_fac->num; j++)
                            zassenhaus_prune_add_factor(Z,
                                  temp_fac->p[j].length - 1, temp_fac->exp[j]);
                        zassenhaus_prune_end_add_factors(Z);

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

        p = (fac->p + 0)->mod.n;
            
        if (r == 1 && r <= cutoff)
        {
            fmpz_poly_factor_insert(final_fac, f, exp);
        }
        else if (r > cutoff && use_van_hoeij)
        {
           fmpz_poly_factor_van_hoeij(final_fac, fac, f, exp, p);
        }
        else
        {
            slong a;
            fmpz_t T;
            fmpz_poly_factor_t lifted_fac;

            fmpz_poly_factor_init(lifted_fac);
            fmpz_init(T);

            fmpz_poly_factor_mignotte(T, f);
            /*
               bound adjustment, we multiply true factors (which might be
               monic) by the leading coefficient of f in the implementation
               below
            */
            fmpz_mul(T, T, f->coeffs + f->length - 1);
            fmpz_abs(T, T);
            fmpz_mul_ui(T, T, 2);
            fmpz_add_ui(T, T, 1);
            a = fmpz_clog_ui(T, p);

            fmpz_poly_hensel_lift_once(lifted_fac, f, fac, a);

            #if TRACE_ZASSENHAUS == 1
            flint_printf("|p = %wd, a = %wd\n", p, a);
            flint_printf("|Pre hensel lift factorisation (nmod_poly):\n");
            nmod_poly_factor_print(fac);
            flint_printf("|Post hensel lift factorisation (fmpz_poly):\n");
            fmpz_poly_factor_print(lifted_fac);
            #endif

            fmpz_set_ui(T, p);
            fmpz_pow_ui(T, T, a);
            fmpz_poly_factor_zassenhaus_recombination_with_prune(
                                          final_fac, lifted_fac, f, T, exp, Z);

            fmpz_poly_factor_clear(lifted_fac);
            fmpz_clear(T);
        }

        nmod_poly_factor_clear(fac);
        zassenhaus_prune_clear(Z);
    }
}

void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    const slong lenG = G->length;
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
        slong j, k;
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
        {
            _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j],
                                                sq_fr_fac->p + j, WORD_MAX, 0);
        }

        fmpz_poly_factor_clear(sq_fr_fac);
    }
    fmpz_poly_clear(g);
}

#undef TRACE_ZASSENHAUS

