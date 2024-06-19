/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_asinh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_asinh_series(g, h, hlen, n, ctx));
}

void
nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_asinh_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_asinh_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}

void
_nmod_poly_asin_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_asin_series(g, h, hlen, n, ctx));
}

void
nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_asin_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_asin_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}

void
_nmod_poly_atanh_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_atanh_series(g, h, hlen, n, ctx));
}

void
nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_atanh_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_atanh_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}

void
_nmod_poly_atan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_atan_series(g, h, hlen, n, ctx));
}

void
nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_atan_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_atan_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}

void
_nmod_poly_cosh_series(nn_ptr f, nn_srcptr h, slong n, nmod_t mod)
{
    nn_ptr g = _nmod_vec_init(n);
    _nmod_poly_exp_expinv_series(f, g, h, n, n, mod);
    _nmod_vec_add(f, f, g, n, mod);
    _nmod_vec_scalar_mul_nmod(f, f, n, n_invmod(UWORD(2), mod.n), mod);
    _nmod_vec_clear(g);
}

void
nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    nn_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    slong h_len;

    h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_cosh_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        if (n > 0)
            nmod_poly_set_coeff_ui(g, 0, UWORD(1));
        return;
    }

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && h_len >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        g_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(g, n);
        g_coeffs = g->coeffs;
    }

    _nmod_poly_cosh_series(g_coeffs, h_coeffs, n, h->mod);

    if (h == g && h_len >= n)
    {
        nmod_poly_swap(g, t1);
        nmod_poly_clear(t1);
    }

    g->length = n;

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(g);
}

void
_nmod_poly_cos_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod)
{
    nn_ptr t, u;

    t = _nmod_vec_init(n);
    u = _nmod_vec_init(n);

    /* cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2) */
    _nmod_vec_scalar_mul_nmod(u, h, n, n_invmod(UWORD(2), mod.n), mod);
    _nmod_poly_tan_series(t, u, n, n, mod);
    _nmod_poly_mullow(u, t, n, t, n, n, mod);
    _nmod_vec_neg(t, u, n, mod);
    t[0] = u[0] = UWORD(1);
    _nmod_poly_div_series(g, t, n, u, n, n, mod);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    nn_ptr h_coeffs;
    slong h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_cos_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        if (n > 0)
            nmod_poly_set_coeff_ui(g, 0, UWORD(1));
        return;
    }

    nmod_poly_fit_length(g, n);

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    _nmod_poly_cos_series(g->coeffs, h_coeffs, n, h->mod);

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    g->length = n;
	_nmod_poly_normalise(g);
}

void
_nmod_poly_sinh_series(nn_ptr f, nn_srcptr h, slong n, nmod_t mod)
{
    nn_ptr g = _nmod_vec_init(n);
    _nmod_poly_exp_expinv_series(f, g, h, n, n, mod);
    _nmod_vec_sub(f, f, g, n, mod);
    _nmod_vec_scalar_mul_nmod(f, f, n, n_invmod(UWORD(2), mod.n), mod);
    _nmod_vec_clear(g);
}

void
nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    nn_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    slong h_len;

    h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_sinh_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        return;
    }

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && h_len >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        g_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(g, n);
        g_coeffs = g->coeffs;
    }

    _nmod_poly_sinh_series(g_coeffs, h_coeffs, n, h->mod);

    if (h == g && h_len >= n)
    {
        nmod_poly_swap(g, t1);
        nmod_poly_clear(t1);
    }

    g->length = n;

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(g);
}

void
_nmod_poly_sin_series(nn_ptr g, nn_srcptr h, slong n, nmod_t mod)
{
    nn_ptr t, u;

    t = _nmod_vec_init(n);
    u = _nmod_vec_init(n);

    /* sin(x) = 2*tan(x/2)/(1+tan(x/2)^2) */
    _nmod_vec_scalar_mul_nmod(u, h, n, n_invmod(2, mod.n), mod);
    _nmod_poly_tan_series(t, u, n, n, mod);
    _nmod_poly_mullow(u, t, n, t, n, n, mod); u[0] = UWORD(1);
    _nmod_poly_div_series(g, t, n, u, n, n, mod);
    _nmod_vec_add(g, g, g, n, mod);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    nn_ptr h_coeffs;
    slong h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_sin_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    _nmod_poly_sin_series(g->coeffs, h_coeffs, n, h->mod);

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    g->length = n;
	_nmod_poly_normalise(g);
}

void
_nmod_poly_tanh_series(nn_ptr f, nn_srcptr h, slong n, nmod_t mod)
{
    nn_ptr t, u;

    t = _nmod_vec_init(n);
    u = _nmod_vec_init(n);

    _nmod_vec_add(t, h, h, n, mod);
    _nmod_poly_exp_series(u, t, n, n, mod);
    _nmod_vec_set(t, u, n);
    t[0] = UWORD(0);
    u[0] = UWORD(2);
    _nmod_poly_div_series(f, t, n, u, n, n, mod);

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void
nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    nn_ptr h_coeffs;
    slong h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_tanh_series). Constant term != 0.\n");
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    _nmod_poly_tanh_series(g->coeffs, h_coeffs, n, h->mod);

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    g->length = n;
	_nmod_poly_normalise(g);
}

void
_nmod_poly_tan_series(nn_ptr g, nn_srcptr h, slong hlen, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_tan_series(g, h, hlen, n, ctx));
}

void
nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, slong n)
{
    slong hlen = h->length;

    if (hlen > 0 && h->coeffs[0] != UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_tan_series). Constant term != 0.\n");
    }

    if (hlen <= 1 || n <= 1)
    {
        nmod_poly_zero(g);
        return;
    }

    nmod_poly_fit_length(g, n);
    _nmod_poly_tan_series(g->coeffs, h->coeffs, hlen, n, h->mod);
    _nmod_poly_set_length(g, n);
	_nmod_poly_normalise(g);
}
