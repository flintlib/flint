/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_HYPGEOM_IMPL_H
#define ACB_HYPGEOM_IMPL_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

slong acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec, slong n_max);

void acb_hypgeom_gamma_stirling_inner(acb_t s, const acb_t z, slong N, slong prec);

int acb_hypgeom_u_asymp_determine_region(const mag_t r, const mag_t zlo, const acb_t z);

void acb_hypgeom_mag_chi(mag_t chi, ulong n);

void _acb_hypgeom_legendre_q_single(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);
void _acb_hypgeom_legendre_q_double(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);

int _mag_gt_norm_ui(const mag_t a, const mag_t b, const mag_t c, ulong n);

#ifdef __cplusplus
}
#endif

#endif /* ACB_HYPGEOM_IMPL_H */
