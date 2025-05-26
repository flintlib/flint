/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* todo: make it possible to compute a superset of the result, implement the
 * algorithm of Gerhard, Giesbrecht, Storjohann, Zima (2003) */

#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"

int
gr_poly_dispersion_resultant(fmpz_t disp, gr_vec_t disp_set,
			     const gr_poly_t f, const gr_poly_t g,
			     gr_ctx_t ctx)
{
    gr_ctx_t Pol_a, ZZ;
    gr_poly_t a, fa, ga, res;
    gr_vec_t roots, mult;

    truth_t is_zero = truth_and(gr_poly_is_zero(f, ctx),
				gr_poly_is_zero(g, ctx));
    if (is_zero == T_TRUE)
	return GR_DOMAIN;
    else if (is_zero == T_UNKNOWN)
	return GR_UNABLE;

    gr_ctx_init_gr_poly(Pol_a, ctx);
    gr_ctx_init_fmpz(ZZ);

    gr_poly_init(a, ctx);
    gr_poly_init(fa, Pol_a);
    gr_poly_init(ga, Pol_a);
    gr_poly_init(res, ZZ);

    gr_vec_init(roots, 0, ZZ);
    gr_vec_init(mult, 0, ZZ);

    int status = GR_SUCCESS;

    status |= gr_poly_set_gr_poly_other(fa, f, ctx, Pol_a);
    status |= gr_poly_set_gr_poly_other(ga, g, ctx, Pol_a);

    status |= gr_poly_gen(a, ctx);
    status |= gr_poly_taylor_shift(fa, fa, a, Pol_a);

    /* todo: composed sum */
    status |= gr_poly_resultant(res, fa, ga, Pol_a);

    status |= gr_poly_roots_other(roots, mult, res, ctx, 0, ZZ);

    if (disp_set != NULL)
	gr_vec_set_length(disp_set, 0, ZZ);

    for (slong i = 0; i < roots->length; i++)
    {
	gr_srcptr rt = gr_vec_entry_srcptr(roots, i, ZZ);
	if (fmpz_sgn(rt) >= 0) {
	    if (disp_set != NULL)
		status |= gr_vec_append(disp_set, rt, ZZ);
	    if (disp != NULL && fmpz_cmp(rt, disp) >= 0)
		fmpz_swap(disp, (fmpz *) rt);
	}
    }

    if (disp_set != NULL)
	_fmpz_vec_sort(disp_set->entries, disp_set->length);

    gr_vec_clear(mult, ZZ);
    gr_vec_clear(roots, ZZ);

    gr_poly_clear(res, ctx);
    gr_poly_clear(ga, Pol_a);
    gr_poly_clear(fa, Pol_a);
    gr_poly_clear(a, ctx);

    gr_ctx_clear(ZZ);
    gr_ctx_clear(Pol_a);

    return status;
}

int
gr_poly_dispersion_from_factors(fmpz_t disp, gr_vec_t disp_set,
			   const gr_vec_t ffac, const gr_vec_t gfac,
			   gr_ctx_t ctx)
{
    gr_ctx_t ZZ, pctx;
    gr_ptr a, b, _alpha;
    fmpz_t alpha;
    gr_poly_t pshift;

    if (gr_ctx_is_unique_factorization_domain(ctx) != T_TRUE)
	return GR_UNABLE;

    /* todo: generalize to positive characteristic? (cf. Gerhard et al.) */
    if (gr_ctx_is_finite_characteristic(ctx) != T_FALSE)
	return GR_UNABLE;

    gr_ctx_init_fmpz(ZZ);

    if (disp_set != NULL)
	gr_vec_set_length(disp_set, 0, ZZ);

    /* We assume that f and g are nonzero, so no factors means constant */
    if (ffac->length == 0 || gfac->length == 0)
	return GR_SUCCESS;

    gr_ctx_init_gr_poly(pctx, ctx);
    GR_TMP_INIT3(a, b, _alpha, ctx);
    fmpz_init(alpha);
    gr_poly_init(pshift, ctx);

    int status = GR_SUCCESS;

    if (ffac == gfac)
    {
	if (disp_set != NULL)
	    status |= gr_vec_append(disp_set, alpha, ZZ);
	if (disp != NULL)
	    fmpz_zero(disp);
    }

    for (slong i = 0; i < ffac->length; i++)
    {
	const gr_poly_struct *p = gr_vec_entry_srcptr(ffac, i, pctx);

	slong j0 = (ffac == gfac) ? i + 1 : 0;
	for (slong j = j0; j < gfac->length; j++)
	{
	    const gr_poly_struct *q = gr_vec_entry_srcptr(gfac, j, pctx);

	    int status1 = GR_SUCCESS;

	    status1 |= gr_poly_leading_taylor_shift(_alpha, p, q, ctx);
	    if (status1 == GR_DOMAIN)
		continue;
	    status1 |= gr_get_fmpz(alpha, _alpha, ctx);
	    if (status1 == GR_DOMAIN)
		continue;
	    status |= status1;

	    if (status != GR_SUCCESS)
		goto cleanup;

	    if (ffac == gfac)
		fmpz_abs(alpha, alpha);
	    else if (fmpz_sgn(alpha) < 0)
		continue;

	    if (disp_set != NULL
		? gr_vec_contains(disp_set, alpha, ZZ) == T_TRUE
		: disp != NULL && fmpz_cmp(alpha, disp) <= 0)
		    continue;

	    if (p->length > 2)  /* GR_SUCCESS => degree well-defined */
	    {
		status |= gr_poly_taylor_shift(pshift, p, _alpha, ctx);
		int eq = gr_poly_equal(pshift, q, ctx);
		if (eq == T_FALSE)
		    continue;
		else if (eq == T_UNKNOWN)
		    status = GR_UNABLE;
	    }

	    if (status != GR_SUCCESS)
		goto cleanup;

	    if (disp_set != NULL)
		status |= gr_vec_append(disp_set, alpha, ZZ);
	    if (disp != NULL && fmpz_cmp(alpha, disp) > 0)
		fmpz_swap(disp, alpha);
	}
    }

    if (disp_set != NULL)
	_fmpz_vec_sort(disp_set->entries, disp_set->length);

cleanup:

    gr_poly_clear(pshift, ctx);
    fmpz_clear(alpha);
    GR_TMP_CLEAR3(a, b, _alpha, ctx);
    gr_ctx_clear(pctx);
    gr_ctx_clear(ZZ);

    return status;
}

int
gr_poly_dispersion_factor(fmpz_t disp, gr_vec_t disp_set,
			  const gr_poly_t f, const gr_poly_t g,
			  gr_ctx_t ctx)
{
    gr_ctx_t pctx, ZZ;
    gr_poly_t fsqf, gsqf;
    gr_vec_t ffac, gfac, ignored;
    gr_ptr c;

    switch (truth_and(gr_poly_is_zero(f, ctx), gr_poly_is_zero(g, ctx)))
    {
	case T_TRUE:
	    return GR_DOMAIN;
	case T_UNKNOWN:
	    return GR_UNABLE;
	case T_FALSE:
	    ;
    }

    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(pctx, ctx);
    gr_poly_init(fsqf, ctx);
    gr_poly_init(gsqf, ctx);
    gr_vec_init(ffac, 0, pctx);
    gr_vec_init(gfac, 0, pctx);
    gr_vec_init(ignored, 0, ZZ);
    GR_TMP_INIT(c, pctx);

    if (gr_poly_squarefree_part(fsqf, f, ctx) != GR_SUCCESS)
	status |= gr_poly_set(fsqf, f, ctx);
    status |= gr_factor(c, ffac, ignored, f, 0, pctx);

    if (f == g)
    {
	status |= gr_poly_dispersion_from_factors(disp, disp_set, ffac, ffac, ctx);
    }
    else
    {
	if (gr_poly_squarefree_part(gsqf, g, ctx) != GR_SUCCESS)
	    status |= gr_poly_set(gsqf, g, ctx);
	status |= gr_factor(c, gfac, ignored, g, 0, pctx);
	status |= gr_poly_dispersion_from_factors(disp, disp_set, ffac, gfac, ctx);
    }

    GR_TMP_CLEAR(c, pctx);
    gr_vec_clear(gfac, pctx);
    gr_vec_clear(ffac, pctx);
    gr_poly_clear(fsqf, ctx);
    gr_poly_clear(gsqf, ctx);
    gr_vec_clear(ignored, ZZ);
    gr_ctx_clear(pctx);
    gr_ctx_clear(ZZ);

    return status;
}

int
gr_poly_dispersion(fmpz_t disp, gr_vec_t disp_set,
		   const gr_poly_t f, const gr_poly_t g,
		   gr_ctx_t ctx)
{
    return gr_poly_dispersion_factor(disp, disp_set, f, g, ctx);
}
