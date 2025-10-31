/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_IMPL_H
#define GR_IMPL_H

#include "nmod_types.h"
#include "fmpz_mod_types.h"
#include "gr_types.h"
#include "qqbar.h"

int _gr_fmpz_poly_factor(fmpz_poly_t c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx);

void _gr_fmpz_mpoly_ctx_clear(gr_ctx_t ctx);
int _gr_fmpz_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);
slong _gr_fmpz_mpoly_ctx_ngens(slong * ngens, gr_ctx_t ctx);
int _gr_fmpz_mpoly_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx);

int gr_ctx_init_fq_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);
int gr_ctx_init_fq_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_nmod_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);
int gr_ctx_init_fq_nmod_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_zech_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);
int gr_ctx_init_fq_zech_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);

int _gr_arf_cmpabs(int * res, const arf_t x, const arf_t y, const gr_ctx_t ctx);
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);

int _gr_arb_cmpabs(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx);

int _gr_acb_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);

int _gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);
int _gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int _gr_gr_poly_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx);

void qqbar_set_fmpzi(qqbar_t res, const fmpzi_t x);

#endif
