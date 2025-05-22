/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "gr_poly.h"

int
gr_poly_add_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0) {
        if (gr_is_zero(c, ctx) == T_TRUE) {
            return gr_poly_zero(res, ctx);
        } else {
            gr_poly_fit_length(res, 1, ctx);
            _gr_poly_set_length(res, 1, ctx);
            status = gr_set(res->coeffs, c, ctx);
            _gr_poly_normalise(res, ctx);
            return status;
        }
    }

    status = GR_SUCCESS;

    if (res != poly) {
        status |= gr_poly_set(res, poly, ctx);
    }

    if (gr_is_zero(c, ctx) != T_TRUE) {
        gr_ptr constant_coeff = gr_poly_coeff_ptr(res, 0, ctx);
        status |= gr_add(constant_coeff, constant_coeff, c, ctx);
        if (len == 1 && gr_is_zero(constant_coeff, ctx) == T_TRUE) {
            _gr_poly_set_length(res, 0, ctx);
        }
    }

    return status;
}

int
gr_poly_add_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0) {
        if (c == 0) {
            return gr_poly_zero(res, ctx);
        } else {
            gr_poly_fit_length(res, 1, ctx);
            _gr_poly_set_length(res, 1, ctx);
            status = gr_set_ui(res->coeffs, c, ctx);
            _gr_poly_normalise(res, ctx);
            return status;
        }
    }

    status = GR_SUCCESS;

    if (res != poly) {
        status |= gr_poly_set(res, poly, ctx);
    }

    if (c != 0) {
        gr_ptr constant_coeff = gr_poly_coeff_ptr(res, 0, ctx);
        status |= gr_add_ui(constant_coeff, constant_coeff, c, ctx);
        if (len == 1 && gr_is_zero(constant_coeff, ctx) == T_TRUE) {
            _gr_poly_set_length(res, 0, ctx);
        }
    }

    return status;
}

int
gr_poly_add_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0) {
        if (c == 0) {
            return gr_poly_zero(res, ctx);
        } else {
            gr_poly_fit_length(res, 1, ctx);
            _gr_poly_set_length(res, 1, ctx);
            status = gr_set_si(res->coeffs, c, ctx);
            _gr_poly_normalise(res, ctx);
            return status;
        }
    }

    status = GR_SUCCESS;

    if (res != poly) {
        status |= gr_poly_set(res, poly, ctx);
    }

    if (c != 0) {
        gr_ptr constant_coeff = gr_poly_coeff_ptr(res, 0, ctx);
        status |= gr_add_si(constant_coeff, constant_coeff, c, ctx);
        if (len == 1 && gr_is_zero(constant_coeff, ctx) == T_TRUE) {
            _gr_poly_set_length(res, 0, ctx);
        }
    }

    return status;
}

int
gr_poly_add_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0) {
        if (fmpz_is_zero(c)) {
            return gr_poly_zero(res, ctx);
        } else {
            gr_poly_fit_length(res, 1, ctx);
            _gr_poly_set_length(res, 1, ctx);
            status = gr_set_fmpz(res->coeffs, c, ctx);
            _gr_poly_normalise(res, ctx);
            return status;
        }
    }

    status = GR_SUCCESS;

    if (res != poly) {
        status |= gr_poly_set(res, poly, ctx);
    }

    if (!fmpz_is_zero(c)) {
        gr_ptr constant_coeff = gr_poly_coeff_ptr(res, 0, ctx);
        status |= gr_add_fmpz(constant_coeff, constant_coeff, c, ctx);
        if (len == 1 && gr_is_zero(constant_coeff, ctx) == T_TRUE) {
            _gr_poly_set_length(res, 0, ctx);
        }
    }

    return status;
}

int
gr_poly_add_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0) {
        if (fmpq_is_zero(c)) {
            return gr_poly_zero(res, ctx);
        } else {
            gr_poly_fit_length(res, 1, ctx);
            _gr_poly_set_length(res, 1, ctx);
            status = gr_set_fmpq(res->coeffs, c, ctx);
            _gr_poly_normalise(res, ctx);
            return status;
        }
    }

    status = GR_SUCCESS;

    if (res != poly) {
        status |= gr_poly_set(res, poly, ctx);
    }

    if (!fmpq_is_zero(c)) {
        gr_ptr constant_coeff = gr_poly_coeff_ptr(res, 0, ctx);
        status |= gr_add_fmpq(constant_coeff, constant_coeff, c, ctx);
        if (len == 1 && gr_is_zero(constant_coeff, ctx) == T_TRUE) {
            _gr_poly_set_length(res, 0, ctx);
        }
    }

    return status;
}
