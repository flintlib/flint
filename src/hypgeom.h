/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef HYPGEOM_H
#define HYPGEOM_H

#include "fmpz_types.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz_poly_t A;
    fmpz_poly_t B;
    fmpz_poly_t P;
    fmpz_poly_t Q;

    /* precomputation data */
    int have_precomputed;
    slong r;
    slong boundC;
    slong boundD;
    slong boundK;
    mag_t MK;
}
hypgeom_struct;

typedef hypgeom_struct hypgeom_t[1];

void hypgeom_init(hypgeom_t hyp);

void hypgeom_clear(hypgeom_t hyp);

void hypgeom_precompute(hypgeom_t hyp);

slong hypgeom_estimate_terms(const mag_t z, int r, slong prec);

slong hypgeom_bound(mag_t error, int r,
    slong C, slong D, slong K, const mag_t TK, const mag_t z, slong prec);

void arb_hypgeom_sum(arb_t P, arb_t Q, const hypgeom_t hyp, slong n, slong prec);

void arb_hypgeom_infsum(arb_t P, arb_t Q, hypgeom_t hyp, slong target_prec, slong prec);

#ifdef __cplusplus
}
#endif

#endif

