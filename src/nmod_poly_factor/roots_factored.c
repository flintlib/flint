/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "long_extras.h"
#include "ulong_extras.h"

typedef nmod_poly_t n_list_t;

/* lets not generate solutions lists with length longer than LENGTH_LIMIT */
#if FLINT64
#define LENGTH_LIMIT (WORD(1) << 32)
#else
#define LENGTH_LIMIT (WORD(1) << 25)
#endif

/*
    The modulus of b is divisible by the modulus of a.
    Map b via the projection.
*/
static void map_down(nmod_poly_t a, const nmod_poly_t b)
{
    slong i;
    FLINT_ASSERT((b->mod.n % a->mod.n) == 0);
    nmod_poly_fit_length(a, b->length);
    for (i = 0; i < b->length; i++)
        a->coeffs[i] = b->coeffs[i] % a->mod.n;
    a->length = b->length;
    _nmod_poly_normalise(a);
}

/*
    Every lifter needs a Diophantine equation solver:
    Given a, b, c with a >= 0, b > 0, try to push all solutions for x to
        a*x + b*y = c
    with 0 <= x < b.
*/
static int dio_solve(n_list_t v, ulong A, ulong B, ulong C)
{
    int success = 1;
    slong k;
    ulong t, d;
    fmpz_t xstart, xstride, xlength;
    fmpz_t a, b, c;

    fmpz_init(xstart);
    fmpz_init(xstride);
    fmpz_init(xlength);
    fmpz_init_set_ui(a, A);
    fmpz_init_set_ui(b, B);
    fmpz_init_set_ui(c, C);

    fmpz_divides_mod_list(xstart, xstride, xlength, c, a, b);

    k = *xlength;
    if ((!COEFF_IS_MPZ(k)) && (k + v->length < LENGTH_LIMIT))
    {
        nmod_poly_fit_length(v, k + v->length);
        t = fmpz_get_ui(xstart);
        d = fmpz_get_ui(xstride);
        for (; k > 0; k--)
        {
            v->coeffs[v->length] = t;
            v->length++;
            t += d;
        }
    }
    else
    {
        /* too many solutions */
        success = 0;
    }

    fmpz_clear(xstart);
    fmpz_clear(xstride);
    fmpz_clear(xlength);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);

    return success;
}

/* Fill x with the roots of fpk, where f->mod.n is p^k */
static int roots_mod_prime_power(nmod_poly_factor_t x, nmod_poly_t fpk,
                                       ulong p, slong k, int with_multiplicity)
{
    int success = 1;
    slong e1, e2;
    slong i, j, old_length;
    n_list_t x1, x2;
    nmod_poly_t f, dfpk, tf, tr, tq;
    ulong pe1, pe2e1, fprime, mfpe1;

    FLINT_ASSERT(k >= 1);
    FLINT_ASSERT(n_is_probabprime(p));

    nmod_poly_init_mod(tf, fpk->mod);
    nmod_poly_init_mod(tr, fpk->mod);
    nmod_poly_init_mod(tq, fpk->mod);

    nmod_poly_init_mod(dfpk, fpk->mod);
    nmod_poly_derivative(dfpk, fpk);

    nmod_poly_init_mod(x1, fpk->mod);
    nmod_poly_init_mod(x2, fpk->mod);

    nmod_poly_init(f, p);
    map_down(f, fpk);

    /* try to fill x1 with solutions mod p */
    x1->length = 0;
    if (f->length > 0)
    {
        nmod_poly_factor_t r;
        nmod_poly_factor_init(r);

        nmod_poly_roots(r, f, 0);
        nmod_poly_fit_length(x1, r->num);
        for (i = 0; i < r->num; i++)
            x1->coeffs[i] = nmod_neg(r->p[i].coeffs[0], f->mod);
        x1->length = r->num;

        nmod_poly_factor_clear(r);
    }
    else
    {
        if (p >= LENGTH_LIMIT)
        {
            /* too many solutions mod p */
            success = 0;
            goto cleanup;
        }

        nmod_poly_fit_length(x1, p);
        for (i = 0; i < p; i++)
            x1->coeffs[i] = i;
        x1->length = p;
    }

    /* lift roots mod p^e1 to roots mod p^e2 */
    for (e1 = 1; e1 < k; e1 = e2)
    {
        e2 = FLINT_MIN(k, 2*e1);

        pe1 = n_pow(p, e1);
        pe2e1 = n_pow(p, e2 - e1);
        x2->length = 0;
        for (i = 0; i < x1->length; i++)
        {
            mfpe1 = nmod_poly_evaluate_nmod(fpk, x1->coeffs[i]);
            mfpe1 = nmod_neg(mfpe1, fpk->mod);
            FLINT_ASSERT((mfpe1 % pe1) == 0);
            mfpe1 = mfpe1/pe1;

            fprime = nmod_poly_evaluate_nmod(dfpk, x1->coeffs[i]);
            fprime = fprime % pe2e1;

            old_length = x2->length;
            if (!dio_solve(x2, fprime, pe2e1, mfpe1))
            {
                success = 0;
                goto cleanup;
            }

            for (j = old_length; j < x2->length; j++)
            {
                x2->coeffs[j] = x1->coeffs[i] + x2->coeffs[j] * pe1;
            }
        }
        nmod_poly_swap(x1, x2);
    }

    /* fill in roots and multiplicies if wanted */
    nmod_poly_factor_fit_length(x, x1->length);
    for (i = 0; i < x1->length; i++)
    {
        nmod_poly_fit_length(x->p + i, 2);
        x->p[i].mod = fpk->mod;          /* bummer */
        x->p[i].coeffs[1] = 1;
        FLINT_ASSERT(x1->coeffs[i] < fpk->mod.n);
        x->p[i].coeffs[0] = nmod_neg(x1->coeffs[i], fpk->mod);
        x->p[i].length = 2;
        x->exp[i] = 1;
        if (with_multiplicity)
        {
            if (fpk->length > 0)
            {
                nmod_poly_divrem(tf, tr, fpk, x->p + i);
                FLINT_ASSERT(nmod_poly_is_zero(tr));
                while (nmod_poly_divrem(tq, tr, tf, x->p + i),
                       nmod_poly_is_zero(tr))
                {
                    FLINT_ASSERT(tf->length >= (x->p + i)->length);
                    x->exp[i]++;
                    nmod_poly_swap(tq, tf);
                }
            }
            else
            {
                x->exp[i] = WORD_MAX;
            }
        }
    }
    x->num = x1->length;

cleanup:

    nmod_poly_clear(tf);
    nmod_poly_clear(tr);
    nmod_poly_clear(tq);

    nmod_poly_clear(x1);
    nmod_poly_clear(x2);
    nmod_poly_clear(f);

    nmod_poly_clear(dfpk);

    return success;
}


int nmod_poly_roots_factored(nmod_poly_factor_t x0,
            const nmod_poly_t f, int with_multiplicity, const n_factor_t * fac)
{
    int success = 1;
    slong i, j, k, new_length;
    ulong m;
    nmod_poly_factor_t x1, x2;
    nmod_poly_t fpe;

    if (f->length <= 0)
    {
        flint_throw(FLINT_ERROR, "Exception in nmod_poly_roots_factored: "
                                                  "input polynomial is zero.");
        return 0;
    }

    nmod_poly_init(fpe, fac->p[0]);

    m = 1;

    nmod_poly_factor_init(x1);
    nmod_poly_factor_init(x2);

    i = 0;
    nmod_poly_init(fpe, n_pow(fac->p[i], fac->exp[i]));
    map_down(fpe, f);
    if (!roots_mod_prime_power(x0, fpe, fac->p[i], fac->exp[i],
                                                            with_multiplicity))
    {
        goto almost_failed;
    }

    for (i = 1; x0->num > 0 && i < fac->num; i++)
    {
        m *= fpe->mod.n;

        nmod_init(&fpe->mod, n_pow(fac->p[i], fac->exp[i]));
        map_down(fpe, f);
        if (!roots_mod_prime_power(x1, fpe, fac->p[i], fac->exp[i],
                                                            with_multiplicity))
        {
            goto almost_failed;
        }

        if (z_mul_checked(&new_length, x0->num, x1->num) ||
            new_length >= LENGTH_LIMIT)
        {
            goto almost_failed;
        }

        /* combine roots with CRT, multiplicities with FLINT_MIN */
        x2->num = 0;
        nmod_poly_factor_fit_length(x2, new_length);
        for (j = 0; j < x0->num; j++)
        for (k = 0; k < x1->num; k++)
        {
            nmod_poly_struct * r = x2->p + x2->num;
            nmod_poly_fit_length(r, 2);
            r->mod = f->mod;        /* bummer */
            r->coeffs[1] = 1;
            FLINT_ASSERT(x1->p[k].length == 2);
            FLINT_ASSERT(x0->p[j].length == 2);
            r->coeffs[0] = n_CRT(x1->p[k].coeffs[0], fpe->mod.n,
                                 x0->p[j].coeffs[0], m);
            r->length = 2;

            FLINT_ASSERT(x0->exp[j] >= 1);
            FLINT_ASSERT(x1->exp[k] >= 1);
            x2->exp[x2->num] = FLINT_MIN(x0->exp[j], x1->exp[k]);
            x2->num++;
        }
        nmod_poly_factor_swap(x0, x2);
    }

cleanup:

    nmod_poly_factor_clear(x1);
    nmod_poly_factor_clear(x2);

    nmod_poly_clear(fpe);

    return success;

almost_failed:

    /* if any prime power is lacking roots, we can still succeed */

    x0->num = 0;

    for (i++; i < fac->num; i++)
    {
        nmod_init(&fpe->mod, n_pow(fac->p[i], fac->exp[i]));
        map_down(fpe, f);
        if (roots_mod_prime_power(x1, fpe, fac->p[i], fac->exp[i], 0) &&
            x1->num == 0)
        {
            goto cleanup;
        }
    }

    success = 0;
    goto cleanup;
}
