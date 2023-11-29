/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2016 Vincent Delecroix
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void
fmpz_poly_randtest_irreducible1(fmpz_poly_t p, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    slong i;
    fmpz_t c;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t q;

    len = 1 + n_randint(state, len);

    fmpz_init(c);

    if (bits == 1)
        fmpz_set_ui(c, 2);
    else
        fmpz_randprime(c, state, bits, 0);
    fmpz_mod_ctx_init(ctx, c);
    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_randtest_irreducible(q, state, len, ctx);

    fmpz_mod_poly_get_fmpz_poly(p, q, ctx);

    /* After lifting, the coefficients belong to {0, ..., c-1}. We now  */
    /* randomly subtract c so that some of them become negative.        */
    for (i = 0; i < p->length; i++)
    {
        if (n_randint(state, 3) == 0 && !(bits == 1 && fmpz_is_zero(p->coeffs + i)))
            fmpz_sub(p->coeffs + i, p->coeffs + i, c);
    }

    fmpz_poly_content(c, p);
    fmpz_poly_scalar_divexact_fmpz(p, p, c);

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(c);
}

void
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
fmpz_poly_randtest_irreducible(fmpz_poly_t pol, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    if (n_randint(state, 2))
        fmpz_poly_randtest_irreducible1(pol, state, len, bits);
    else
        fmpz_poly_randtest_irreducible2(pol, state, len, bits);
}

void
fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state,
                   slong len, flint_bitcnt_t bits)
{
    fmpz_poly_fit_length(f, len);
    _fmpz_vec_randtest(f->coeffs, state, len, bits);
    _fmpz_poly_set_length(f, len);
    _fmpz_poly_normalise(f);
}

void
fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state,
                            slong len, flint_bitcnt_t bits)
{
    fmpz_poly_fit_length(f, len);
    _fmpz_vec_randtest_unsigned(f->coeffs, state, len, bits);
    _fmpz_poly_set_length(f, len);
    _fmpz_poly_normalise(f);
}

void
fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state,
                            slong len, flint_bitcnt_t bits)
{
    if ((bits == 0) || (len == 0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_randtest_not_zero). bits or len is zero.\n");
    }

    fmpz_poly_randtest(f, state, len, bits);
    if (fmpz_poly_is_zero(f))
        fmpz_poly_set_ui(f, 1);
}
