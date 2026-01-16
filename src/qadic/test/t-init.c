/*
    Copyright (C) 2026 Alexey Orlov

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "qadic.h"
#include "nmod_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/*
originally it was possible to create qadic context only with conway polynomial
  or, using qadic_ctx_init, with a random irreducible polynomial if it fails

we have two cubic irreducible polynomials over F_2
x^3 + x + 1 (conway polynomial)
x^3 + x^2 + 1
we create 4 qadic contexts:
  1. with conway polynomial, using _qadic_ctx_init_conway_ui
  2. using qadic_ctx_init, which will also use conway polynomial
  3. with nmod_poly_t modulus f = x^3 + x^2 + 1
  4. with fmpz_mod_poly_t modulus f = x^3 + x^2 + 1
for all 4 we compute g^3 + g^2 + 1 and g^3 + g + 1 for a generator g (lift of a generator of the residue field)
*/

static inline void compute_f(qadic_t out, qadic_ctx_t ctx)
{
    qadic_t z, z0, z2, z3, tmp;

    qadic_init(z); qadic_init(z0); qadic_init(z2); qadic_init(z3); qadic_init(tmp);

    qadic_gen(z, ctx); qadic_set_ui(z0, 1, ctx); qadic_mul(z2, z, z, ctx); qadic_mul(z3, z2, z, ctx);
    qadic_add(tmp, z0, z2, ctx);
    qadic_add(out, tmp, z3, ctx);

    qadic_clear(z); qadic_clear(z0); qadic_clear(z2); qadic_clear(z3); qadic_clear(tmp);
}

static inline void compute_f_conway(qadic_t out, qadic_ctx_t ctx)
{
    qadic_t z, z0, z2, z3, tmp;

    qadic_init(z); qadic_init(z0); qadic_init(z2); qadic_init(z3); qadic_init(tmp);

    qadic_gen(z, ctx); qadic_set_ui(z0, 1, ctx); qadic_mul(z2, z, z, ctx); qadic_mul(z3, z2, z, ctx);
    qadic_add(tmp, z0, z, ctx);
    qadic_add(out, tmp, z3, ctx);

    qadic_clear(z); qadic_clear(z0); qadic_clear(z2); qadic_clear(z3); qadic_clear(tmp);
}

static inline void error_exit(const char* prefix, qadic_t val, qadic_ctx_t ctx)
{
    flint_printf("FAIL\n\n");
    flint_printf("%s", prefix), qadic_print_pretty(val, ctx), flint_printf("\n");
    fflush(stdout);
    flint_abort();
}

TEST_FUNCTION_START(qadic_init, state)
{
    fmpz_mod_ctx_t GF;
    nmod_poly_t f_nmod; fmpz_mod_poly_t f;
    qadic_ctx_t ctx_conway; qadic_ctx_t ctx; qadic_ctx_t ctx_f_nmod; qadic_ctx_t ctx_f;
    qadic_t fg_val;
    fmpz_t p;

    fmpz_init(p); fmpz_set_ui(p, 2);

    nmod_poly_init(f_nmod, 2);
    nmod_poly_set_coeff_ui(f_nmod, 0, 1);
    nmod_poly_set_coeff_ui(f_nmod, 2, 1);
    nmod_poly_set_coeff_ui(f_nmod, 3, 1);

    fmpz_mod_ctx_init_ui(GF, 2);
    fmpz_mod_poly_init(f, GF);
    fmpz_mod_poly_set_nmod_poly(f, f_nmod);

    _qadic_ctx_init_conway_ui(ctx_conway, 2, 3, 0, 1, "x", PADIC_SERIES);
    qadic_ctx_init(ctx, p, 3, 0, 1, "x", PADIC_SERIES);
    qadic_ctx_init_modulus_nmod(ctx_f_nmod, 2, f_nmod, 0, 1, "x", PADIC_SERIES);
    qadic_ctx_init_modulus(ctx_f, p, f, 0, 1, "x", PADIC_SERIES);

    qadic_init(fg_val);

    compute_f(fg_val, ctx_conway);
    if (qadic_is_zero(fg_val))
        error_exit("conway: g^3 + g^2 + 1 = ", fg_val, ctx_conway);

    compute_f_conway(fg_val, ctx_conway);
    if (!qadic_is_zero(fg_val))
        error_exit("conway: g^3 + g + 1 = ", fg_val, ctx_conway);

    compute_f(fg_val, ctx);
    if (qadic_is_zero(fg_val))
        error_exit("default: g^3 + g^2 + 1 = ", fg_val, ctx);

    compute_f_conway(fg_val, ctx);
    if (!qadic_is_zero(fg_val))
        error_exit("default: g^3 + g + 1 = ", fg_val, ctx);

    compute_f(fg_val, ctx_f);
    if (!qadic_is_zero(fg_val))
        error_exit("f_fmpz: g^3 + g^2 + 1 = ", fg_val, ctx_f);

    compute_f_conway(fg_val, ctx_f);
    if (qadic_is_zero(fg_val))
        error_exit("f_fmpz: g^3 + g + 1 = ", fg_val, ctx_f);

    compute_f(fg_val, ctx_f_nmod);
    if (!qadic_is_zero(fg_val))
        error_exit("f_nmod: g^3 + g^2 + 1 = ", fg_val, ctx_f_nmod);

    compute_f_conway(fg_val, ctx_f_nmod);
    if (qadic_is_zero(fg_val))
        error_exit("f_nmod: g^3 + g + 1 = ", fg_val, ctx_f_nmod);

    qadic_clear(fg_val);
    qadic_ctx_clear(ctx_conway);
    qadic_ctx_clear(ctx);
    qadic_ctx_clear(ctx_f_nmod);
    qadic_ctx_clear(ctx_f);
    nmod_poly_clear(f_nmod);
    fmpz_mod_poly_clear(f, GF);
    fmpz_mod_ctx_clear(GF);
    fmpz_clear(p);

    TEST_FUNCTION_END(state);
}

