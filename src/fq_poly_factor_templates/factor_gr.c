/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "templates.h"

/*
    Wrappers around the gr_poly implementations. A polynomial of this
    module has the same layout as a gr_poly_struct, so polynomials can be
    cast; the factors are moved (not copied) out of the gr_poly_vec.

    Note: the gr context constructed here refers to (rather than owns)
    the given context, so it must not be cleared.
*/

static void
TEMPLATE(T, poly_factor_set_gr) (TEMPLATE(T, poly_factor_t) res,
    gr_poly_vec_t fac, const fmpz_vec_t exp, gr_ctx_t gr_ctx,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, num = fac->length;

    TEMPLATE(T, poly_factor_fit_length) (res, res->num + num, ctx);

    for (i = 0; i < num; i++)
    {
        TEMPLATE(T, poly_struct) * p = res->poly + res->num + i;

        TEMPLATE(T, poly_clear) (p, ctx);
        *((gr_poly_struct *) p) = *(fac->entries + i);
        gr_poly_init(fac->entries + i, gr_ctx);

        res->exp[res->num + i] = fmpz_get_si(exp->entries + i);
    }

    res->num += num;
}

/* Also used by factor_berlekamp.c, factor_cantor_zassenhaus.c and
   factor_kaltofen_shoup.c, which this file precedes in the
   instantiations. */
static void
TEMPLATE(T, poly_factor_gr_algorithm) (TEMPLATE(T, poly_factor_t) res,
    TEMPLATE(T, t) leading_coeff, const TEMPLATE(T, poly_t) input,
    int algorithm, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    gr_ptr lc;

    res->num = 0;

    if (input->length == 0)
    {
        if (leading_coeff != NULL)
            TEMPLATE(T, zero) (leading_coeff, ctx);
        return;
    }

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);
    lc = gr_heap_init(gr_ctx);

    GR_MUST_SUCCEED(_gr_poly_factor_finite_field(lc, fac, exp,
        (const gr_poly_struct *) input, algorithm, gr_ctx));

    if (leading_coeff != NULL)
        TEMPLATE(T, set) (leading_coeff, lc, ctx);

    TEMPLATE(T, poly_factor_set_gr) (res, fac, exp, gr_ctx, ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
    gr_heap_clear(lc, gr_ctx);
}

void
TEMPLATE(T, poly_factor) (TEMPLATE(T, poly_factor_t) res,
    TEMPLATE(T, t) leading_coeff, const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_factor_gr_algorithm) (res, leading_coeff, input,
        GR_POLY_FACTOR_ALGORITHM_DEFAULT, ctx);
}

void
TEMPLATE(T, poly_factor_with_berlekamp) (TEMPLATE(T, poly_factor_t) res,
    TEMPLATE(T, t) leading_coeff, const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_factor_gr_algorithm) (res, leading_coeff, input,
        GR_POLY_FACTOR_ALGORITHM_BERLEKAMP, ctx);
}

void
TEMPLATE(T, poly_factor_with_cantor_zassenhaus) (TEMPLATE(T, poly_factor_t) res,
    TEMPLATE(T, t) leading_coeff, const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_factor_gr_algorithm) (res, leading_coeff, input,
        GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS, ctx);
}

void
TEMPLATE(T, poly_factor_with_kaltofen_shoup) (TEMPLATE(T, poly_factor_t) res,
    TEMPLATE(T, t) leading_coeff, const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_factor_gr_algorithm) (res, leading_coeff, input,
        GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP, ctx);
}

void
TEMPLATE(T, poly_factor_squarefree) (TEMPLATE(T, poly_factor_t) res,
    const TEMPLATE(T, poly_t) f, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    gr_ptr lc;

    res->num = 0;

    if (f->length <= 1)
        return;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);
    lc = gr_heap_init(gr_ctx);

    GR_MUST_SUCCEED(gr_poly_factor_squarefree(lc, fac, exp,
        (const gr_poly_struct *) f, gr_ctx));

    TEMPLATE(T, poly_factor_set_gr) (res, fac, exp, gr_ctx, ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
    gr_heap_clear(lc, gr_ctx);
}

void
TEMPLATE(T, poly_factor_distinct_deg) (TEMPLATE(T, poly_factor_t) res,
    const TEMPLATE(T, poly_t) poly, slong * const * degs,
    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t degrees;
    slong i, num;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(degrees, 0);

    GR_MUST_SUCCEED(gr_poly_factor_distinct_deg(fac, degrees,
        (const gr_poly_struct *) poly, gr_ctx));

    num = fac->length;
    TEMPLATE(T, poly_factor_fit_length) (res, res->num + num, ctx);

    for (i = 0; i < num; i++)
    {
        TEMPLATE(T, poly_struct) * p = res->poly + res->num + i;

        TEMPLATE(T, poly_clear) (p, ctx);
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
TEMPLATE(T, poly_factor_equal_deg) (TEMPLATE(T, poly_factor_t) factors,
    const TEMPLATE(T, poly_t) pol, slong d, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    slong i;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_poly_vec_init(fac, 0, gr_ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(gr_poly_factor_equal_deg(fac,
        (const gr_poly_struct *) pol, d, gr_ctx));

    for (i = 0; i < fac->length; i++)
        fmpz_vec_append_ui(exp, 1);

    TEMPLATE(T, poly_factor_set_gr) (factors, fac, exp, gr_ctx, ctx);

    gr_poly_vec_clear(fac, gr_ctx);
    fmpz_vec_clear(exp);
}

/* Note: unlike gr_poly_is_irreducible, this module considers constants
   (and the zero polynomial) irreducible. */
static int
TEMPLATE(T, poly_irreducible_gr) (const TEMPLATE(T, poly_t) f, int ddf,
    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    truth_t res;

    if (f->length <= 2)
        return 1;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    res = ddf ? gr_poly_is_irreducible_ddf((const gr_poly_struct *) f, gr_ctx)
              : gr_poly_is_irreducible((const gr_poly_struct *) f, gr_ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "poly_is_irreducible: unable to decide\n");

    return (res == T_TRUE);
}

int
TEMPLATE(T, poly_is_irreducible) (const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, poly_irreducible_gr) (f, 0, ctx);
}

int
TEMPLATE(T, poly_is_irreducible_ddf) (const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, poly_irreducible_gr) (f, 1, ctx);
}

void
TEMPLATE(T, poly_roots) (TEMPLATE(T, poly_factor_t) r,
    const TEMPLATE(T, poly_t) f, int with_multiplicity,
    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_vec_t roots;
    fmpz_vec_t mult;
    slong i, num, sz;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    /* fq contexts already know they are fields */
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_vec_init(roots, 0, gr_ctx);
    fmpz_vec_init(mult, 0);

    GR_MUST_SUCCEED(gr_poly_roots_finite_field(roots, mult,
        (const gr_poly_struct *) f, 0, gr_ctx));

    num = roots->length;
    sz = gr_ctx->sizeof_elem;

    r->num = 0;
    TEMPLATE(T, poly_factor_fit_length) (r, num, ctx);

    for (i = 0; i < num; i++)
    {
        TEMPLATE(T, poly_struct) * p = r->poly + i;
        TEMPLATE(T, t) c;

        TEMPLATE(T, init) (c, ctx);
        TEMPLATE(T, neg) (c, (TEMPLATE(T, struct) *) GR_ENTRY(roots->entries, i, sz), ctx);
        TEMPLATE(T, poly_zero) (p, ctx);
        TEMPLATE(T, poly_set_coeff) (p, 0, c, ctx);
        TEMPLATE(T, one) (c, ctx);
        TEMPLATE(T, poly_set_coeff) (p, 1, c, ctx);
        TEMPLATE(T, clear) (c, ctx);

        r->exp[i] = with_multiplicity ? fmpz_get_si(mult->entries + i) : 1;
    }

    r->num = num;

    gr_vec_clear(roots, gr_ctx);
    fmpz_vec_clear(mult);
}


/* Finds one linear factor of a polynomial that splits into linear
   factors (which need not be distinct). */
void
TEMPLATE(T, poly_factor_split_single) (TEMPLATE(T, poly_t) linfactor,
    const TEMPLATE(T, poly_t) input, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_vec_t roots;
    fmpz_vec_t mult;
    TEMPLATE(T, t) c;

    if (input->length == 2)
    {
        TEMPLATE(T, poly_set) (linfactor, input, ctx);
        return;
    }

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    gr_vec_init(roots, 0, gr_ctx);
    fmpz_vec_init(mult, 0);

    GR_MUST_SUCCEED(gr_poly_roots_finite_field(roots, mult,
        (const gr_poly_struct *) input, 0, gr_ctx));

    if (roots->length == 0)
        flint_throw(FLINT_ERROR, "poly_factor_split_single: no linear factor\n");

    TEMPLATE(T, init) (c, ctx);
    TEMPLATE(T, neg) (c, (TEMPLATE(T, struct) *) roots->entries, ctx);
    TEMPLATE(T, poly_zero) (linfactor, ctx);
    TEMPLATE(T, poly_set_coeff) (linfactor, 0, c, ctx);
    TEMPLATE(T, one) (c, ctx);
    TEMPLATE(T, poly_set_coeff) (linfactor, 1, c, ctx);
    TEMPLATE(T, clear) (c, ctx);

    gr_vec_clear(roots, gr_ctx);
    fmpz_vec_clear(mult);
}

int
TEMPLATE(T, poly_is_irreducible_ben_or) (const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    truth_t res;

    if (f->length <= 2)
        return 1;

    TEMPLATE3(_gr_ctx_init, T, from_ref) (gr_ctx, ctx);
    GR_IGNORE(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    res = gr_poly_is_irreducible_ben_or((const gr_poly_struct *) f, gr_ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "poly_is_irreducible_ben_or: unable to decide\n");

    return (res == T_TRUE);
}

#endif
