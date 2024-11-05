/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_IMPL_H
#define CA_IMPL_H

#include "ca_types.h"
#include "fexpr.h"

#ifdef __cplusplus
extern "C" {
#endif

int ca_ext_can_evaluate_qqbar(const ca_ext_t x, ca_ctx_t ctx);

void ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx);
void _ca_all_extensions(ca_ext_ptr ** extensions, slong * length, const ca_t x, ca_ctx_t ctx);

void
fexpr_set_nf_elem(fexpr_t res, const nf_elem_t a, const nf_t nf, const fexpr_t var);

void
_ca_get_fexpr_given_ext(fexpr_t res, const ca_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void
_ca_ext_get_fexpr_given_ext(fexpr_t res, const ca_ext_t x, ulong flags,
        ca_ext_ptr * ext, slong num_ext, const fexpr_struct * ext_vars, ca_ctx_t ctx);

void
_ca_default_variables(fexpr_ptr ext_vars, slong num_ext);

void
_nf_elem_get_fmpz_poly_den_shallow(fmpz_poly_t pol, fmpz_t den, const nf_elem_t a, const nf_t nf);

void
ca_set_ext(ca_t res, ca_ext_srcptr ext, ca_ctx_t ctx);

void
ca_clear_unchecked(ca_t x, ca_ctx_t ctx);

ca_field_ptr ca_ctx_get_field_qqbar(ca_ctx_t ctx, const qqbar_t x);

int
_ca_cmp(const ca_t x, const ca_t y, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* CA_IMPL_H */
