/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

/*
    Wrappers around the gr_poly implementations. An fmpz_mod_poly_struct
    has the same layout as a gr_poly_struct, so polynomials can be cast;
    the factors are moved (not copied) out of the gr_poly_vec.
*/
static void
_fmpz_mod_poly_factor_set_gr(fmpz_mod_poly_factor_t res, gr_poly_vec_t fac,
    const fmpz_vec_t exp, const fmpz_mod_ctx_t ctx, gr_ctx_t gr_ctx)
{
    slong i, num = fac->length;

    fmpz_mod_poly_factor_fit_length(res, res->num + num, ctx);

    for (i = 0; i < num; i++)
    {
        fmpz_mod_poly_struct * p = res->poly + res->num + i;

        fmpz_mod_poly_clear(p, ctx);
        *((gr_poly_struct *) p) = *(fac->entries + i);
        gr_poly_init(fac->entries + i, gr_ctx);

        res->exp[res->num + i] = fmpz_get_si(exp->entries + i);
    }

    res->num += num;
}

static void
_fmpz_mod_poly_factor_gr(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f,
    int algorithm, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    fmpz_t lc;

    res->num = 0;

    if (f->length <= 1)
        return;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);
    fmpz_init(lc);

    GR_MUST_SUCCEED(_gr_poly_factor_finite_field(lc, fac, exp,
        (const gr_poly_struct *) f, algorithm, gr_ctx));

    _fmpz_mod_poly_factor_set_gr(res, fac, exp, ctx, gr_ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
    fmpz_clear(lc);
}

void
fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_t f,
    const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_DEFAULT, ctx);
}

void
fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS, ctx);
}

void
fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_BERLEKAMP, ctx);
}

void
fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP, ctx);
}

void
fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    fmpz_t lc;

    res->num = 0;

    if (f->length <= 1)
        return;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);
    fmpz_init(lc);

    GR_MUST_SUCCEED(gr_poly_factor_squarefree(lc, fac, exp,
        (const gr_poly_struct *) f, gr_ctx));

    _fmpz_mod_poly_factor_set_gr(res, fac, exp, ctx, gr_ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
    fmpz_clear(lc);
}

void
fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t poly, slong * const * degs, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t degrees;
    slong i, num;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(degrees, 0);

    GR_MUST_SUCCEED(gr_poly_factor_distinct_deg(fac, degrees,
        (const gr_poly_struct *) poly, gr_ctx));

    num = fac->length;
    fmpz_mod_poly_factor_fit_length(res, res->num + num, ctx);

    for (i = 0; i < num; i++)
    {
        fmpz_mod_poly_struct * p = res->poly + res->num + i;

        fmpz_mod_poly_clear(p, ctx);
        *((gr_poly_struct *) p) = *(fac->entries + i);
        gr_poly_init(fac->entries + i, gr_ctx);

        res->exp[res->num + i] = 1;
        (*degs)[i] = fmpz_get_si(degrees->entries + i);
    }

    res->num += num;

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(degrees);
}

void
fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors,
    const fmpz_mod_poly_t f, slong d, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    slong i;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(gr_poly_factor_equal_deg(fac, (const gr_poly_struct *) f, d, gr_ctx));

    for (i = 0; i < fac->length; i++)
        fmpz_vec_append_ui(exp, 1);

    _fmpz_mod_poly_factor_set_gr(factors, fac, exp, ctx, gr_ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
}

/* Note: unlike gr_poly_is_irreducible, this module considers constants
   (and the zero polynomial) irreducible. */
static int
_fmpz_mod_poly_irreducible_gr(const fmpz_mod_poly_t f, int ddf, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    truth_t res;

    if (f->length <= 2)
        return 1;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    res = ddf ? gr_poly_is_irreducible_ddf((const gr_poly_struct *) f, gr_ctx)
              : gr_poly_is_irreducible((const gr_poly_struct *) f, gr_ctx);


    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "fmpz_mod_poly_is_irreducible: unable to decide\n");

    return (res == T_TRUE);
}

int
fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_irreducible_gr(f, 0, ctx);
}

int
fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_irreducible_gr(f, 1, ctx);
}

void
fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f,
    int with_mult, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_vec_t roots;
    fmpz_vec_t mult;
    slong i, num;

    /* the context borrows ctx and must not be cleared */
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_vec_init(roots, 0, gr_ctx);
    fmpz_vec_init(mult, 0);

    GR_MUST_SUCCEED(gr_poly_roots_finite_field(roots, mult,
        (const gr_poly_struct *) f, 0, gr_ctx));

    num = roots->length;
    r->num = 0;
    fmpz_mod_poly_factor_fit_length(r, num, ctx);

    for (i = 0; i < num; i++)
    {
        fmpz_mod_poly_struct * p = r->poly + i;
        fmpz_t c;

        fmpz_init(c);
        fmpz_mod_neg(c, ((const fmpz *) roots->entries) + i, ctx);
        fmpz_mod_poly_zero(p, ctx);
        fmpz_mod_poly_set_coeff_fmpz(p, 0, c, ctx);
        fmpz_mod_poly_set_coeff_ui(p, 1, 1, ctx);
        fmpz_clear(c);

        r->exp[i] = with_mult ? fmpz_get_si(mult->entries + i) : 1;
    }

    r->num = num;

    gr_vec_clear(roots, gr_ctx);
    fmpz_vec_clear(mult);
}

/* The distinct degree factorization uses several threads by itself when
   they are available (see gr_poly_factor_distinct_deg). */
void
fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t poly, slong * const * degs, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_factor_distinct_deg(res, poly, degs, ctx);
}
