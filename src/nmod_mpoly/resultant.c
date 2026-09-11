/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "nmod_mpoly.h"
#include "nmod_mpoly_factor.h"

/* The dense bivariate algorithm works on the whole bounding box of the
   inputs, so it is only worth converting to when the inputs fill enough of
   theirs. Both of them must, and the box must stay within a size that can be
   held at once. */
#define BIVARIATE_DENSITY 8
#define BIVARIATE_MAX_AREA (WORD(1) << 24)

/* res_var(A, B) through the dense bivariate algorithm, when the inputs are
   bivariate and dense enough for it to pay. Returns 0 when they are not, or
   when the algorithm declines. */
static int
_nmod_mpoly_resultant_bivariate(nmod_mpoly_t R, const nmod_mpoly_t A,
            const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong other = 1 - var;
    slong degAy, degAx, degBy, degBx, lenA, lenB, npoints;
    n_bpoly_t Ab, Bb, Rb;
    n_poly_t res;
    int success;

    if (ctx->minfo->nvars != 2 || A->length == 0 || B->length == 0)
        return 0;

    degAy = nmod_mpoly_degree_si(A, var, ctx);
    degAx = nmod_mpoly_degree_si(A, other, ctx);
    degBy = nmod_mpoly_degree_si(B, var, ctx);
    degBx = nmod_mpoly_degree_si(B, other, ctx);

    if (degAy < 1 || degBy < 1)
        return 0;

    lenA = degAy + 1;
    lenB = degBy + 1;

    if ((double) lenA * (degAx + 1) > (double) BIVARIATE_MAX_AREA ||
        (double) lenB * (degBx + 1) > (double) BIVARIATE_MAX_AREA)
        return 0;

    if ((double) A->length * BIVARIATE_DENSITY < (double) lenA * (degAx + 1) ||
        (double) B->length * BIVARIATE_DENSITY < (double) lenB * (degBx + 1))
        return 0;

    npoints = degBy * degAx + degAy * degBx + 1;

    if (!n_bpoly_mod_resultant_cutoff(FLINT_MAX(lenA, lenB),
                                      FLINT_MIN(lenA, lenB), npoints))
        return 0;

    n_bpoly_init(Ab);
    n_bpoly_init(Bb);
    n_bpoly_init(Rb);
    n_poly_init(res);

    /* var becomes the outer variable, the one the resultant eliminates */
    nmod_mpoly_get_bpoly(Ab, A, var, other, ctx);
    nmod_mpoly_get_bpoly(Bb, B, var, other, ctx);

    success = n_bpoly_mod_resultant(res, Ab, Bb, ctx->mod);

    if (success)
    {
        /* the resultant is a polynomial in the other variable alone */
        n_bpoly_fit_length(Rb, 1);
        n_poly_swap(Rb->coeffs + 0, res);
        Rb->length = 1;
        n_bpoly_normalise(Rb);

        nmod_mpoly_set_bpoly(R, A->bits, Rb, var, other, ctx);
    }

    n_poly_clear(res);
    n_bpoly_clear(Ab);
    n_bpoly_clear(Bb);
    n_bpoly_clear(Rb);

    return success;
}

int nmod_mpoly_resultant(nmod_mpoly_t R, const nmod_mpoly_t A,
           const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpoly_univar_t Ax, Bx;

    if (_nmod_mpoly_resultant_bivariate(R, A, B, var, ctx))
        return 1;

    nmod_mpoly_univar_init(Ax, ctx);
    nmod_mpoly_univar_init(Bx, ctx);

    nmod_mpoly_to_univar(Ax, A, var, ctx);
    nmod_mpoly_to_univar(Bx, B, var, ctx);

    success = nmod_mpoly_univar_resultant(R, Ax, Bx, ctx);

    nmod_mpoly_univar_clear(Ax, ctx);
    nmod_mpoly_univar_clear(Bx, ctx);

    return success;
}
