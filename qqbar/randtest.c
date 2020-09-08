/*
    Copyright (C) 2016 Vincent Delecroix
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

#if __FLINT_RELEASE >= 20700

static void
fmpz_poly_randtest_irreducible1(fmpz_poly_t p, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    slong i;
    fmpz_t c;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t q;

    len = 1 + n_randint(state, len);

    fmpz_init(c);

    fmpz_randprime(c, state, bits, 0);
    fmpz_mod_ctx_init(ctx, c);
    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_randtest_irreducible(q, state, len, ctx);

    fmpz_mod_poly_get_fmpz_poly(p, q, ctx);

    /* After lifting, the coefficients belong to {0, ..., c-1}. We now  */
    /* randomly subtract c so that some of them become negative.        */
    for (i = 0; i < p->length; i++)
    {
        if (n_randint(state, 3) == 0)
            fmpz_sub(
                fmpz_poly_get_coeff_ptr(p, i),
                fmpz_poly_get_coeff_ptr(p, i),
                c);
    }

    fmpz_poly_content(c, p);
    fmpz_poly_scalar_divexact_fmpz(p, p, c);

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(c);
}

#else

static void
fmpz_poly_randtest_irreducible1(fmpz_poly_t p, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    slong i;
    fmpz_t c;
    fmpz_mod_poly_t q;

    len = 1 + n_randint(state, len);

    fmpz_init(c);

    fmpz_randprime(c, state, bits, 0);
    fmpz_mod_poly_init(q, c);
    fmpz_mod_poly_randtest_irreducible(q, state, len);

    fmpz_mod_poly_get_fmpz_poly(p, q);

    /* After lifting, the coefficients belong to {0, ..., c-1}. We now  */
    /* randomly subtract c so that some of them become negative.        */
    for (i = 0; i < p->length; i++)
    {
        if (n_randint(state, 3) == 0)
            fmpz_sub(
                fmpz_poly_get_coeff_ptr(p, i),
                fmpz_poly_get_coeff_ptr(p, i),
                c);
    }

    fmpz_poly_content(c, p);
    fmpz_poly_scalar_divexact_fmpz(p, p, c);

    fmpz_mod_poly_clear(q);
    fmpz_clear(c);
}

#endif

static void
fmpz_poly_randtest_irreducible2(fmpz_poly_t pol, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    while (1)
    {
        slong i;
        fmpz_poly_factor_t fac;

        do {
            fmpz_poly_randtest(pol, state, len, bits);
        } while (fmpz_poly_degree(pol) < 1);

        fmpz_poly_factor_init(fac);
        fmpz_poly_factor(fac, pol);

        i = n_randint(state, fac->num);

        if (FLINT_ABS(fmpz_poly_max_bits(fac->p + i)) <= bits)
        {
            fmpz_poly_set(pol, fac->p + i);
            fmpz_poly_factor_clear(fac);
            break;
        }

        fmpz_poly_factor_clear(fac);
    }
}

void
_qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits, int real)
{
    fmpz_poly_t pol;
    slong prec, i, rdeg, r1, r2;
    acb_ptr roots;

    deg = FLINT_MAX(deg, 1);
    bits = FLINT_MAX(bits, 1);

    if (deg == 1 || n_randint(state, 4) == 0)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_randtest(t, state, bits);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    fmpz_poly_init(pol);

    do {
        if (n_randint(state, 2))
            fmpz_poly_randtest_irreducible1(pol, state, deg + 1, bits);
        else
            fmpz_poly_randtest_irreducible2(pol, state, deg + 1, bits);

        rdeg = fmpz_poly_degree(pol);
        r1 = rdeg;
        r2 = 0;
        if (real)
            fmpz_poly_signature(&r1, &r2, pol);
    }
    while (rdeg < 1 || (real == 1 && r1 < 1) || (real == 2 && r2 < 1));

    if (fmpz_sgn(pol->coeffs + rdeg) < 0)
        fmpz_poly_neg(pol, pol);

    roots = _acb_vec_init(rdeg);
    if (real == 0)
        i = n_randint(state, rdeg);
    else if (real == 1)
        i = n_randint(state, r1);
    else
        i = r1 + n_randint(state, 2 * r2);

    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        arb_fmpz_poly_complex_roots(roots, pol, 0, prec);

        if (_qqbar_validate_uniqueness(roots + i, pol, roots + i, 2 * prec))
        {
            fmpz_poly_set(QQBAR_POLY(res), pol);
            acb_set(QQBAR_ENCLOSURE(res), roots + i);
            break;
        }
    }

    _acb_vec_clear(roots, rdeg);
    fmpz_poly_clear(pol);
}

void
qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    _qqbar_randtest(res, state, deg, bits, 0);
}

void
qqbar_randtest_real(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    _qqbar_randtest(res, state, deg, bits, 1);
}

void
qqbar_randtest_nonreal(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    if (deg <= 1)
    {
        flint_printf("qqbar_randtest_nonreal: must have deg >= 2\n");
        flint_abort();
    }

    _qqbar_randtest(res, state, deg, bits, 2);
}

