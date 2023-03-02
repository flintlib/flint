/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "fmpz.h"
#include "fmpq.h"
#include "long_extras.h"

typedef fmpz_poly_t fmpz_list_t;

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
static void map_down(fmpz_mod_poly_t a, const fmpz_mod_ctx_t ctxa,
               const fmpz_mod_poly_t b, const fmpz_mod_ctx_t ctxb)
{
    slong i;

    FLINT_ASSERT(fmpz_divisible(fmpz_mod_ctx_modulus(ctxb),
                                fmpz_mod_ctx_modulus(ctxa)));

    fmpz_mod_poly_fit_length(a, b->length, ctxa);
    for (i = 0; i < b->length; i++)
        fmpz_mod(a->coeffs + i, b->coeffs + i, fmpz_mod_ctx_modulus(ctxa));
    a->length = b->length;
    _fmpz_mod_poly_normalise(a);
}


/* Fill x with the roots of fpk, where f->p is p^k */
static int roots_mod_prime_power(fmpz_mod_poly_factor_t x,
                         fmpz_mod_poly_t fpk, slong k, int with_multiplicity,
                         const fmpz_mod_ctx_t ctxp, const fmpz_mod_ctx_t ctxpk)
{
    int success = 1;
    slong i, j, e1, e2;
    fmpz_list_t x1, x2;
    fmpz_mod_poly_t f, dfpk, tf, tr, tq;
    fmpz_t pe1, pe2e1, fprime, mfpe1;
    fmpz_t xstart, xstride, xlength;

    FLINT_ASSERT(k >= 1);
    FLINT_ASSERT(fmpz_is_probabprime(fmpz_mod_ctx_modulus(ctxp)));

    fmpz_init(pe1);
    fmpz_init(pe2e1);
    fmpz_init(fprime);
    fmpz_init(mfpe1);
    fmpz_init(xstart);
    fmpz_init(xstride);
    fmpz_init(xlength);

    fmpz_mod_poly_init(tf, ctxpk);
    fmpz_mod_poly_init(tr, ctxpk);
    fmpz_mod_poly_init(tq, ctxpk);

    fmpz_mod_poly_init(dfpk, ctxpk);
    fmpz_mod_poly_derivative(dfpk, fpk, ctxpk);

    fmpz_poly_init(x1);
    fmpz_poly_init(x2);

    fmpz_mod_poly_init(f, ctxp);

    map_down(f, ctxp, fpk, ctxpk);

    /* try to fill x1 with solution mod p */
    x1->length = 0;
    if (f->length > 0)
    {
        fmpz_mod_poly_factor_t r;
        fmpz_mod_poly_factor_init(r, ctxp);

        fmpz_mod_poly_roots(r, f, 0, ctxp);
        fmpz_poly_fit_length(x1, r->num);
        for (i = 0; i < r->num; i++)
            fmpz_mod_neg(x1->coeffs + i, r->poly[i].coeffs + 0, ctxp);
        x1->length = r->num;
        fmpz_mod_poly_factor_clear(r, ctxp);
    }
    else
    {
        slong p_si;

        if (fmpz_cmp_si(fmpz_mod_ctx_modulus(ctxp), LENGTH_LIMIT) >= 0)
        {
            /* too many solution mod p */
            success = 0;
            goto cleanup;
        }

        p_si = fmpz_get_si(fmpz_mod_ctx_modulus(ctxp));

        fmpz_poly_fit_length(x1, p_si);
        for (i = 0; i < p_si; i++)
            fmpz_set_si(x1->coeffs + i, i);
        x1->length = p_si;
    }

    /* lift roots mod p^e1 to roots mod p^e2 */
    for (e1 = 1; e1 < k; e1 = e2)
    {
        e2 = FLINT_MIN(k, 2*e1);

        fmpz_pow_ui(pe1, fmpz_mod_ctx_modulus(ctxp), e1);
        fmpz_pow_ui(pe2e1, fmpz_mod_ctx_modulus(ctxp), e2 - e1);
        x2->length = 0;
        for (i = 0; i < x1->length; i++)
        {
            fmpz_mod_poly_evaluate_fmpz(mfpe1, fpk, x1->coeffs + i, ctxpk);
            fmpz_neg(mfpe1, mfpe1);
            FLINT_ASSERT(fmpz_divisible(mfpe1, pe1));
            fmpz_divexact(mfpe1, mfpe1, pe1);

            fmpz_mod_poly_evaluate_fmpz(fprime, dfpk, x1->coeffs + i, ctxpk);
            fmpz_mod(fprime, fprime, pe2e1);

            fmpz_divides_mod_list(xstart, xstride, xlength, mfpe1, fprime, pe2e1);

            j = *xlength;
            if ((!COEFF_IS_MPZ(j)) && (j + x2->length < LENGTH_LIMIT))
            {
                fmpz_poly_fit_length(x2, j + x2->length);
                for (; j > 0; j--)
                {
                    FLINT_ASSERT(x2->length < x2->alloc);
                    fmpz_set(x2->coeffs + x2->length, x1->coeffs + i);
                    fmpz_addmul(x2->coeffs + x2->length, xstart, pe1);
                    fmpz_add(xstart, xstart, xstride);
                    x2->length++;
                }
            }
            else
            {
                /* too many solutions */
                success = 0;
                goto cleanup;
            }
        }
        fmpz_poly_swap(x1, x2);
    }

    /* fill in roots and multiplicies if wanted */
    fmpz_mod_poly_factor_fit_length(x, x1->length, ctxpk);
    for (i = 0; i < x1->length; i++)
    {
        fmpz_mod_poly_fit_length(x->poly + i, 2, ctxpk);
        fmpz_one(x->poly[i].coeffs + 1);

        fmpz_mod_neg(x->poly[i].coeffs + 0, x1->coeffs + i, ctxpk);
        _fmpz_mod_poly_set_length(x->poly + i, 2);
        x->exp[i] = 1;
        if (with_multiplicity)
        {
            if (fpk->length > 0)
            {
                fmpz_mod_poly_divrem(tf, tr, fpk, x->poly + i, ctxpk);
                FLINT_ASSERT(fmpz_mod_poly_is_zero(tr, ctxpk));
                while (fmpz_mod_poly_divrem(tq, tr, tf, x->poly + i, ctxpk),
                       fmpz_mod_poly_is_zero(tr, ctxpk))
                {
                    FLINT_ASSERT(tf->length >= (x->poly + i)->length);
                    x->exp[i]++;
                    fmpz_mod_poly_swap(tq, tf, ctxpk);
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

    fmpz_mod_poly_clear(tf, ctxpk);
    fmpz_mod_poly_clear(tr, ctxpk);
    fmpz_mod_poly_clear(tq, ctxpk);

    fmpz_poly_clear(x1);
    fmpz_poly_clear(x2);
    fmpz_mod_poly_clear(f, ctxp);
    fmpz_clear(pe1);
    fmpz_clear(pe2e1);
    fmpz_clear(fprime);
    fmpz_clear(mfpe1);
    fmpz_clear(xstart);
    fmpz_clear(xstride);
    fmpz_clear(xlength);

    fmpz_mod_poly_clear(dfpk, ctxpk);

    return success;
}


int fmpz_mod_poly_roots_factored(fmpz_mod_poly_factor_t x0,
                           const fmpz_mod_poly_t f, int with_mult,
                             const fmpz_factor_t fac, const fmpz_mod_ctx_t ctx)
{
    int success = 1;
    slong i, j, k, new_length;
    fmpz_t m, pe;
    fmpz_mod_poly_factor_t x1, x2;
    fmpz_mod_poly_t fpe;
    fmpz_mod_ctx_t ctxp, ctxpe;

    if (f->length < 1)
    {
        flint_throw(FLINT_ERROR, "Exception in fmpz_mod_poly_roots_factored: "
                                                  "input polynomial is zero.");
        return 0;
    }

    fmpz_init(pe);
    fmpz_mod_poly_init(fpe, ctx);

    fmpz_init_set_ui(m, 1);

    fmpz_mod_poly_factor_init(x1, ctx);
    fmpz_mod_poly_factor_init(x2, ctx);

    i = 0;
    fmpz_mod_ctx_init(ctxp, fac->p + i);
    fmpz_pow_ui(pe, fac->p + i, fac->exp[i]);
    fmpz_mod_ctx_init(ctxpe, pe);
    map_down(fpe, ctxpe, f, ctx);
    if (!roots_mod_prime_power(x0, fpe, fac->exp[i], with_mult, ctxp, ctxpe))
    {
        goto almost_failed;
    }

    for (i = 1; x0->num > 0 && i < fac->num; i++)
    {
        fmpz_mul(m, m, pe);

        fmpz_mod_ctx_set_modulus(ctxp, fac->p + i);
        fmpz_pow_ui(pe, fac->p + i, fac->exp[i]);
        fmpz_mod_ctx_set_modulus(ctxpe, pe);
        map_down(fpe, ctxpe, f, ctx);
        if (!roots_mod_prime_power(x1, fpe, fac->exp[i], with_mult, ctxp, ctxpe))
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
        fmpz_mod_poly_factor_fit_length(x2, new_length, ctx);
        for (j = 0; j < x0->num; j++)
        for (k = 0; k < x1->num; k++)
        {
            fmpz_mod_poly_struct * r = x2->poly + x2->num;
            fmpz_mod_poly_fit_length(r, 2, ctx);
            fmpz_one(r->coeffs + 1);
            FLINT_ASSERT(x1->poly[k].length == 2);
            FLINT_ASSERT(x0->poly[j].length == 2);
            fmpz_CRT(r->coeffs + 0, x1->poly[k].coeffs + 0, pe,
                                    x0->poly[j].coeffs + 0, m, 0);
            _fmpz_mod_poly_set_length(r, 2);

            FLINT_ASSERT(x0->exp[j] >= 1);
            FLINT_ASSERT(x1->exp[k] >= 1);
            x2->exp[x2->num] = FLINT_MIN(x0->exp[j], x1->exp[k]);
            x2->num++;
        }
        fmpz_mod_poly_factor_swap(x0, x2, ctx);
    }

cleanup:

    fmpz_mod_poly_factor_clear(x1, ctx);
    fmpz_mod_poly_factor_clear(x2, ctx);

    fmpz_clear(m);
    fmpz_clear(pe);

    fmpz_mod_poly_clear(fpe, ctx);

    fmpz_mod_ctx_clear(ctxp);
    fmpz_mod_ctx_clear(ctxpe);

    return success;

almost_failed:

    /* if any prime power is lacking roots, we can still succeed */

    x0->num = 0;

    for (i++; i < fac->num; i++)
    {
        fmpz_mod_ctx_set_modulus(ctxp, fac->p + i);
        fmpz_pow_ui(pe, fac->p + i, fac->exp[i]);
        fmpz_mod_ctx_set_modulus(ctxpe, pe);
        map_down(fpe, ctxpe, f, ctx);
        if (roots_mod_prime_power(x1, fpe, fac->exp[i], 0, ctxp, ctxpe) &&
            x1->num == 0)
        {
            goto cleanup;
        }
    }

    success = 0;
    goto cleanup;
}
