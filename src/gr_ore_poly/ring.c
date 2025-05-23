/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

/* Memory management */

void gr_ore_poly_init(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
gr_ore_poly_init2(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
{
    gr_poly_init((gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
    gr_poly_fit_length((gr_poly_struct *) poly, len, GR_ORE_POLY_ELEM_CTX(ctx));
}

void gr_ore_poly_clear(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    gr_poly_clear((gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
}

void
_gr_ore_poly_set_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
{
    _gr_poly_set_length((gr_poly_struct *) poly, len, GR_ORE_POLY_ELEM_CTX(ctx));
}

void
gr_ore_poly_fit_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
{
    gr_poly_fit_length((gr_poly_struct *) poly, len, GR_ORE_POLY_ELEM_CTX(ctx));
}

/* Basic manipulation */

void
_gr_ore_poly_normalise(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    _gr_poly_normalise((gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set((gr_poly_struct *) res, (const gr_poly_struct *) src, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_truncate(gr_ore_poly_t poly, const gr_ore_poly_t src, slong newlen, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_truncate((gr_poly_struct *) poly, (gr_poly_struct *) src, newlen, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_one(gr_ore_poly_t res, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_one((gr_poly_struct *) res, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_neg_one(gr_ore_poly_t res, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_neg_one((gr_poly_struct *) res, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_gen(gr_ore_poly_t res, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_gen((gr_poly_struct *) res, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
_gr_ore_poly_write(gr_stream_t out, gr_srcptr poly, slong n, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_write(out, poly, n, GR_ORE_POLY_CTX(ctx)->var, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_write(gr_stream_t out, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_write(out, poly->coeffs, poly->length, GR_ORE_POLY_CTX(ctx)->var, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
_gr_ore_poly_get_str(char ** res, gr_srcptr f, slong len, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_get_str(res, f, len, GR_ORE_POLY_CTX(ctx)->var, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_get_str(char ** res, const gr_ore_poly_t f, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_get_str(res, f->coeffs, f->length, GR_ORE_POLY_CTX(ctx)->var, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_print(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_ore_poly_write(out, poly, ctx);
}

int
gr_ore_poly_set_str(gr_ore_poly_t res, const char * s, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set_str((gr_poly_struct *) res, s, GR_ORE_POLY_CTX(ctx)->var, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
_gr_ore_poly_set_str(gr_ptr res, const char * s, slong len, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_set_str(res, s, GR_ORE_POLY_CTX(ctx)->var, len, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_randtest(gr_ore_poly_t poly, flint_rand_t state, slong len, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_randtest((gr_poly_struct *) poly, state, len, GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
_gr_ore_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_equal(poly1, len1, poly2, len2, GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_equal(const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_equal((const gr_poly_struct *) poly1, (const gr_poly_struct *) poly2, GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_is_zero(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_is_zero((const gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_is_one(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_is_one((const gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
}

truth_t
gr_ore_poly_is_gen(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_is_gen((const gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set_si(gr_ore_poly_t poly, slong x, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set_si((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set_ui(gr_ore_poly_t poly, ulong x, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set_ui((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set_fmpz(gr_ore_poly_t poly, const fmpz_t x, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set_fmpz((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set_fmpq(gr_ore_poly_t poly, const fmpq_t x, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_set_fmpq((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_set_other(gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)
{
    if (x_ctx == ctx) {
        return gr_poly_set((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx == GR_ORE_POLY_ELEM_CTX(ctx)) {
        return gr_poly_set_scalar((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (GR_ORE_POLY_CTX(ctx)->base_ring->which_ring == GR_CTX_GR_POLY && \
                POLYNOMIAL_CTX(GR_ORE_POLY_CTX(ctx)->base_ring)->base_ring == x_ctx) {
        gr_poly_t tmp;
        tmp->coeffs = (gr_ptr) x;
        tmp->length = 1;
        tmp->alloc = 1;
        return gr_poly_set_scalar((gr_poly_struct *) poly, tmp, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPZ) {
        return gr_poly_set_fmpz((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPQ) {
        return gr_poly_set_fmpq((gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    }
    return GR_UNABLE;
}

/* Arithmetic */

int
gr_ore_poly_neg(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_neg((gr_poly_struct *) res, (const gr_poly_struct *) src, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
_gr_ore_poly_add(gr_ptr res, gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_add(res, poly1, len1, poly2, len2, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_add(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_add((gr_poly_struct *) res, (const gr_poly_struct *) poly1, (const gr_poly_struct *) poly2, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
_gr_ore_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
{
    return _gr_poly_sub((gr_poly_struct *) res, poly1, len1, poly2, len2, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_sub(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_sub((gr_poly_struct *) res, (const gr_poly_struct *) poly1, (const gr_poly_struct *) poly2, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_add_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)
{
    if (x_ctx == ctx) {
        return gr_poly_add((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx == GR_ORE_POLY_ELEM_CTX(ctx)) {
        return gr_poly_add_scalar((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (GR_ORE_POLY_CTX(ctx)->base_ring->which_ring == GR_CTX_GR_POLY && \
                POLYNOMIAL_CTX(GR_ORE_POLY_CTX(ctx)->base_ring)->base_ring == x_ctx) {
        gr_poly_t tmp;
        tmp->coeffs = (gr_ptr) x;
        tmp->length = 1;
        tmp->alloc = 1;
        return gr_poly_add_scalar((gr_poly_struct *) poly, (gr_poly_struct *) poly, tmp, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPZ) {
        return gr_poly_add_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPQ) {
        return gr_poly_add_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    }
    return GR_UNABLE;
}

int
gr_ore_poly_add_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_add_ui((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_add_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_add_si((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_add_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_add_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_add_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_add_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_sub_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)
{
    if (x_ctx == ctx) {
        return gr_poly_sub((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx == GR_ORE_POLY_ELEM_CTX(ctx)) {
        return gr_poly_sub_scalar((gr_poly_struct *) res, (gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (GR_ORE_POLY_CTX(ctx)->base_ring->which_ring == GR_CTX_GR_POLY && \
                POLYNOMIAL_CTX(GR_ORE_POLY_CTX(ctx)->base_ring)->base_ring == x_ctx) {
        gr_poly_t tmp;
        tmp->coeffs = (gr_ptr) x;
        tmp->length = 1;
        tmp->alloc = 1;
        return gr_poly_sub_scalar((gr_poly_struct *) poly, (gr_poly_struct *) poly, tmp, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPZ) {
        return gr_poly_sub_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPQ) {
        return gr_poly_sub_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    }
    return GR_UNABLE;
}

int
gr_ore_poly_sub_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_sub_ui((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_sub_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_sub_si((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_sub_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_sub_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_sub_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_sub_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_other_mul(gr_ore_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    if (x_ctx == ctx) {
        /* TODO: Multiplication of Ore polynomials. */
        return GR_UNABLE;
    } else if (x_ctx == GR_ORE_POLY_ELEM_CTX(ctx)) {
        return gr_poly_scalar_mul((gr_poly_struct *) res, x, (const gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (GR_ORE_POLY_CTX(ctx)->base_ring->which_ring == GR_CTX_GR_POLY && \
            POLYNOMIAL_CTX(GR_ORE_POLY_CTX(ctx)->base_ring)->base_ring == x_ctx) {
        gr_poly_t tmp;
        tmp->coeffs = (gr_ptr) x;
        tmp->length = 1;
        tmp->alloc = 1;
        return gr_poly_scalar_mul((gr_poly_struct *) poly, tmp, (const gr_poly_struct *) poly, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPZ) {
        return gr_poly_mul_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    } else if (x_ctx->which_ring == GR_CTX_FMPQ) {
        return gr_poly_mul_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, x, GR_ORE_POLY_ELEM_CTX(ctx));
    }
    return GR_UNABLE;
}

int
gr_ore_poly_mul_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_mul_ui((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_mul_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_mul_si((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_mul_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_mul_fmpz((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}

int
gr_ore_poly_mul_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
{
    return gr_poly_mul_fmpq((gr_poly_struct *) res, (const gr_poly_struct *) poly, c, GR_ORE_POLY_ELEM_CTX(ctx));
}
