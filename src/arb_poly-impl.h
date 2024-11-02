/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_POLY_IMPL_H
#define ARB_POLY_IMPL_H

#include "arb_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void
_arb_poly_lgamma_series_at_one(arb_ptr u, slong len, slong prec);

void
_arb_poly_gamma_stirling_eval2(arb_ptr res, const arb_t z, slong n, slong num, int diff, slong prec);

void
_arb_poly_gamma_stirling_eval(arb_ptr res, const arb_t z, slong n, slong num, slong prec);

int
_gr_arb_poly_taylor_shift(arb_ptr res, arb_srcptr poly, slong n, const arb_t c, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* ARB_POLY_IMPL_H */
