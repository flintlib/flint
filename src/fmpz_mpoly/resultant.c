/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly_factor.h"

/* The dense bivariate algorithm works on the whole bounding box of the
   inputs, so it is only worth converting to when the inputs fill enough of
   theirs. Both of them must, and the box must stay within a size that can be
   held at once. */
#define BIVARIATE_DENSITY 8
#define BIVARIATE_MAX_AREA (WORD(1) << 22)

/* res_var(A, B) through the dense bivariate multimodular algorithm, when the
   inputs are bivariate and dense enough for it to pay. Returns 0 when they
   are not, or when the algorithm declines. */
static int
_fmpz_mpoly_resultant_bivariate(fmpz_mpoly_t R, const fmpz_mpoly_t A,
            const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong other = 1 - var;
    slong degAy, degAx, degBy, degBx, lenA, lenB, npoints;
    fmpz_bpoly_t Ab, Bb, Rb;
    fmpz_poly_t res;
    int success;

    if (ctx->minfo->nvars != 2 || A->length == 0 || B->length == 0)
        return 0;

    degAy = fmpz_mpoly_degree_si(A, var, ctx);
    degAx = fmpz_mpoly_degree_si(A, other, ctx);
    degBy = fmpz_mpoly_degree_si(B, var, ctx);
    degBx = fmpz_mpoly_degree_si(B, other, ctx);

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

    fmpz_bpoly_init(Ab);
    fmpz_bpoly_init(Bb);
    fmpz_bpoly_init(Rb);
    fmpz_poly_init(res);

    /* var becomes the outer variable, the one the resultant eliminates */
    fmpz_mpoly_get_bpoly(Ab, A, var, other, ctx);
    fmpz_mpoly_get_bpoly(Bb, B, var, other, ctx);

    success = fmpz_bpoly_resultant(res, Ab, Bb, 1);

    if (success)
    {
        /* the resultant is a polynomial in the other variable alone */
        fmpz_bpoly_fit_length(Rb, 1);
        fmpz_poly_swap(Rb->coeffs + 0, res);
        Rb->length = fmpz_poly_is_zero(Rb->coeffs + 0) ? 0 : 1;

        fmpz_mpoly_set_fmpz_bpoly(R, A->bits, Rb, var, other, ctx);
    }

    fmpz_poly_clear(res);
    fmpz_bpoly_clear(Ab);
    fmpz_bpoly_clear(Bb);
    fmpz_bpoly_clear(Rb);

    return success;
}

int fmpz_mpoly_resultant(fmpz_mpoly_t R, const fmpz_mpoly_t A,
                   const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_univar_t Ax, Bx;

    if (_fmpz_mpoly_resultant_bivariate(R, A, B, var, ctx))
        return 1;

    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);

    fmpz_mpoly_to_univar(Ax, A, var, ctx);
    fmpz_mpoly_to_univar(Bx, B, var, ctx);

    success = fmpz_mpoly_univar_resultant(R, Ax, Bx, ctx);

    fmpz_mpoly_univar_clear(Ax, ctx);
    fmpz_mpoly_univar_clear(Bx, ctx);

    return success;
}
