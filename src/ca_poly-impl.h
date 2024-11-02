/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_POLY_IMPL_H
#define CA_POLY_IMPL_H

#include "ca_types.h"
#include "fexpr.h"

#ifdef __cplusplus
extern "C" {
#endif

ca_field_ptr
_ca_vec_same_field2(ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, ca_ctx_t ctx);

void
ca_poly_get_fexpr(fexpr_t res, const ca_poly_t A, ulong flags, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* CA_POLY_IMPL_H */
