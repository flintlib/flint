/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

int _fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char * var)
{
    int ret;
    ulong conway_poly[410]; /* Largest degree in database is 409 */
    nmod_poly_struct tmp;

    ret = _nmod_poly_conway(conway_poly, p, d);

    if (!ret)
        return 0;

    nmod_poly_init(&tmp, p);
    tmp.coeffs = conway_poly;
    tmp.length = d + 1;

    fq_nmod_ctx_init_modulus(ctx, &tmp, var);

    /* No need to clear tmp */

    return 1;
}


void fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char * var)
{
    if (!_fq_nmod_ctx_init_conway_ui(ctx, p, d, var))
        flint_throw(FLINT_ERROR,
                "Exception (fq_nmod_ctx_init_conway_ui).  "
                "The polynomial for (p, d) = (%wu, %wd) is not present in the "
                "database.\n", p, d);

    ctx->is_conway = 1;
}

void fq_nmod_ctx_init_ui(fq_nmod_ctx_t ctx, ulong p, slong d, const char *var)
{
    if (_fq_nmod_ctx_init_conway_ui(ctx, p, d, var))
    {
        ctx->is_conway = 1;
        return;
    }
    else
    {
        flint_rand_t state;
        nmod_poly_t poly;

        ctx->is_conway = 0;

        flint_randinit(state);

        nmod_poly_init2(poly, p, d + 1);
        nmod_poly_randtest_sparse_irreducible(poly, state, d + 1);

        fq_nmod_ctx_init_modulus(ctx, poly, var);

        nmod_poly_clear(poly);
        flint_randclear(state);
    }
}

/* Deprecated functions ******************************************************/

void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_nmod_ctx_init_ui(ctx, fmpz_get_ui(p), d, var); }
int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { return _fq_nmod_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_nmod_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
