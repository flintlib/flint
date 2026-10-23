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

/* The dense bivariate representation holds the whole bounding box of the
   inputs, so these bound what converting to it can cost. They are also what
   keeps the size arithmetic below from overflowing, which is why the
   comparisons are made in double: a box of this many entries is already that
   many words of input before the resultant is even called.

   There is deliberately no density condition. Measured against the sparse
   algorithm, how full the box is barely matters: the subresultant PRS over
   the univariate form pays for the degrees whether or not the coefficients
   are there, so at degree 16 and above the dense algorithm is ahead even at
   a few terms per input. What decides is the degree, which
   n_bpoly_mod_resultant_cutoff already tests. */
/* Lower than over a word-size prime field: the multimodular algorithm holds
   the residues of every prime at once, so its working space is a multiple of
   the number of points rather than of the input alone. */
#define BIVARIATE_MAX_AREA 1073741824.0     /* 2^30 */
#define BIVARIATE_MAX_POINTS 268435456.0    /* 2^28 */

/* A sparse input expands into a box far larger than itself. Below
   BIVARIATE_FREE_AREA entries that does not matter, the box being small
   whatever the input; above it the box may exceed the input by at most
   BIVARIATE_MAX_EXPANSION, so that a few terms of high degree are left to the
   sparse algorithm rather than materialised. */
#define BIVARIATE_FREE_AREA 8388608.0       /* 2^23 */
#define BIVARIATE_MAX_EXPANSION 64.0

/* Two polynomials of only a few terms each have a resultant the subresultant
   PRS finds almost at once, and at low degree the dense algorithm cannot make
   back the cost of the whole box. Measured over Z at degree 8, four terms
   each: three to a hundred times slower, at every coefficient size from 32 to
   10000 bits. From BIVARIATE_MIN_AREA in the degrees up it wins again, the
   PRS by then paying for the degrees whether or not the coefficients are
   there. */
#define BIVARIATE_MIN_TERMS 5
#define BIVARIATE_MIN_AREA 256

/* res_var(A, B) through the dense bivariate multimodular algorithm, when the
   inputs are bivariate and dense enough for it to pay. Returns 0 when they
   are not, or when the algorithm declines. */
static int
_fmpz_mpoly_resultant_bivariate(fmpz_mpoly_t R, const fmpz_mpoly_t A,
            const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong other = 1 - var;
    slong degAy, degAx, degBy, degBx, lenA, lenB, npoints;
    double areaA, areaB;
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

    areaA = (double) lenA * (degAx + 1);
    areaB = (double) lenB * (degBx + 1);

    if (areaA > BIVARIATE_MAX_AREA || areaB > BIVARIATE_MAX_AREA)
        return 0;

    if ((double) degBy * degAx + (double) degAy * degBx + 1.0
            > BIVARIATE_MAX_POINTS)
        return 0;

    if (areaA > BIVARIATE_FREE_AREA &&
            areaA > BIVARIATE_MAX_EXPANSION * A->length)
        return 0;

    if (areaB > BIVARIATE_FREE_AREA &&
            areaB > BIVARIATE_MAX_EXPANSION * B->length)
        return 0;

    if (FLINT_MIN(A->length, B->length) < BIVARIATE_MIN_TERMS &&
            lenA * lenB < BIVARIATE_MIN_AREA)
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
