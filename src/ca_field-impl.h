/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_FIELD_IMPL_H
#define CA_FIELD_IMPL_H

#include "ca_types.h"

#ifdef __cplusplus
extern "C" {
#endif

slong ca_depth(const ca_t x, ca_ctx_t ctx);

ca_field_ptr ca_field_cache_lookup_qqbar(ca_field_cache_t cache, const qqbar_t x, ca_ctx_t ctx);

void _ca_field_ideal_insert_clear_mpoly(ca_field_t K, fmpz_mpoly_t poly, fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);

void fmpz_mpoly_set_coeff_si_x(fmpz_mpoly_t poly, slong c, slong x_var, slong x_exp, const fmpz_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* CA_FIELD_IMPL_H */
