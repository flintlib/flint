/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_POLY_IMPL_H
#define ACB_POLY_IMPL_H

#include "acb_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int polylog_is_real(const acb_t s, const acb_t z);

int
_gr_acb_poly_taylor_shift(acb_ptr res, acb_srcptr poly, slong n, const acb_t c, gr_ctx_t ctx);

void _acb_poly_gamma_stirling_eval(acb_ptr res, const acb_t z, slong n, slong num, slong prec);
void _acb_poly_gamma_stirling_eval2(acb_ptr res, const acb_t z, slong n, slong num, int diff, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* ACB_POLY_IMPL_H */
