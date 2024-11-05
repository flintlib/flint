/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_IMPL_H
#define GR_IMPL_H

#include "fmpz_types.h"
#include "qqbar.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* todo: move */
void qqbar_set_fmpzi(qqbar_t res, const fmpzi_t x);

int _gr_arf_cmpabs(int * res, const arf_t x, const arf_t y, const gr_ctx_t ctx);

int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);

int _gr_arb_cmpabs(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx);

/* hidden feature: also works with arb ctx */
int _gr_acb_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);

int _gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);
int _gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int gr_ctx_init_fq_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);
int gr_ctx_init_fq_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);

int gr_ctx_init_fq_nmod_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_nmod_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);

int gr_ctx_init_fq_zech_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_zech_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);

void _gr_fmpz_mpoly_ctx_clear(gr_ctx_t ctx);
int _gr_fmpz_mpoly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s);

int gr_test_pow_ui_base_scalar_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_ui_base_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_ui_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags);

int gr_test_pow_ui_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_fmpz_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags);

#ifdef __cplusplus
}
#endif

#endif /* GR_IMPL_H */
