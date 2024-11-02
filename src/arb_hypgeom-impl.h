/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_HYPGEOM_IMPL_H
#define ARB_HYPGEOM_IMPL_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int _arf_increment_fast(arf_t x, slong prec);
void _arb_increment_fast(arb_t x, slong prec);

double arf_get_d_log2_abs_approx_clamped(const arf_t x);

void
acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec);

void
arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec);

void
arb_hypgeom_gamma_stirling_inner(arb_t s, const arb_t z, slong N, slong prec);

int
arb_hypgeom_gamma_exact(arb_t res, const arb_t x, int reciprocal, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* ARB_HYPGEOM_IMPL_H */
