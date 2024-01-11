/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

extern int flint_conway_polynomials [];

int _fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char *var)
{
    unsigned int position;

    if (p > 109987)
    {
        return 0;
    }

    for (position = 0; flint_conway_polynomials[position] != 0; position += 3+flint_conway_polynomials[position+1])
    {
        /* Different prime? */
        if (p != flint_conway_polynomials[position])
            continue;

        /* Same degree? */
        if (d == flint_conway_polynomials[position+1])
        {
            nmod_poly_t mod;
            slong i;

            nmod_poly_init(mod, p);

            /* Copy the polynomial */
            for (i = 0; i < d; i++)
            {
                int coeff = flint_conway_polynomials[position+2+i];
                nmod_poly_set_coeff_ui(mod, i, coeff);
            }

            nmod_poly_set_coeff_ui(mod, d, 1);

            fq_nmod_ctx_init_modulus(ctx, mod, var);

            nmod_poly_clear(mod);
            return 1;
        }
    }

    return 0;
}


void fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char *var)
{
    int result;

    if (p > 109987)
    {
        flint_throw(FLINT_ERROR, "Exception (fq_nmod_ctx_init_conway_ui).  Conway polynomials are only available for primes up to 109987.\n");
    }

    result = _fq_nmod_ctx_init_conway_ui(ctx, p, d, var);
    if (!result)
        flint_throw(FLINT_ERROR, "Exception (fq_nmod_ctx_init_conway_ui).  The polynomial for \n"
                "(p,d) = (%wu,%wd) is not present in the database.\n", p, d);

    ctx->is_conway = 1;
}

void fq_nmod_ctx_init_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char *var)
{
    flint_rand_t state;
    nmod_poly_t poly;

    if (_fq_nmod_ctx_init_conway_ui(ctx, p, d, var))
    {
        ctx->is_conway = 1;
	return;
    } else
        ctx->is_conway = 0;

    flint_randinit(state);

    nmod_poly_init2(poly, p, d + 1);
    nmod_poly_randtest_sparse_irreducible(poly, state, d + 1);

    fq_nmod_ctx_init_modulus(ctx, poly, var);

    nmod_poly_clear(poly);
    flint_randclear(state);
}

/* Deprecated functions ******************************************************/

void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_nmod_ctx_init_ui(ctx, fmpz_get_ui(p), d, var); }
int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { return _fq_nmod_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_nmod_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
