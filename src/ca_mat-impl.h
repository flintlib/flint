/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_MAT_IMPL_H
#define CA_MAT_IMPL_H

#include "ca_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void
_ca_set_nf_fmpz_poly_den(ca_t res, const fmpz_poly_t poly, const fmpz_t den, ca_field_t K, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* CA_MAT_IMPL_H */
