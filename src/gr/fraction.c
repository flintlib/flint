/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_generic.h"

typedef struct
{
    gr_ctx_struct * domain;
    int flags;
} _gr_fraction_ctx_struct;

typedef gr_ctx_struct gr_fraction_ctx_struct;
typedef gr_fraction_ctx_struct gr_fraction_ctx_t[1];


#define GR_FRACTION_CTX(ring_ctx) ((_gr_fraction_ctx_struct *)((ring_ctx)))
#define GR_FRACTION_DOMAIN_CTX(ring_ctx) (GR_FRACTION_CTX(ring_ctx)->domain)
#define GR_FRACTION_FLAGS(ring_ctx) (GR_FRACTION_CTX(ring_ctx)->flags)

#define DOMAIN(ctx) GR_FRACTION_DOMAIN_CTX(ctx)
#define NUMER(x, ctx) (x)
#define DENOM(x, ctx) GR_ENTRY((x), 1, GR_FRACTION_DOMAIN_CTX(ctx)->sizeof_elem)






int
_gr_fraction_ctx_write(gr_stream_t out, gr_fraction_ctx_t ctx)
{
    gr_stream_write(out, "Fraction field over ");
    gr_ctx_write(out, GR_FRACTION_DOMAIN_CTX(ctx));
    return GR_SUCCESS;
}

void
_gr_fraction_ctx_clear(gr_ctx_t ctx)
{
}

/* Normally the domain is assumed to be an integral domain, but it could
   be something inexact like polynomials with floating-point coefficients. */
truth_t _gr_fraction_ctx_is_certainly_field_else_unknown(gr_fraction_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(GR_FRACTION_DOMAIN_CTX(ctx)) == T_TRUE ? T_TRUE : T_UNKNOWN;
}

truth_t _gr_fraction_ctx_is_threadsafe(gr_ctx_t ctx) { return gr_ctx_is_threadsafe(GR_FRACTION_DOMAIN_CTX(ctx)); }
truth_t _gr_fraction_ctx_is_finite(gr_ctx_t ctx) { return gr_ctx_is_finite(GR_FRACTION_DOMAIN_CTX(ctx)); }
truth_t _gr_fraction_ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(GR_FRACTION_DOMAIN_CTX(ctx)); }
truth_t _gr_fraction_ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(GR_FRACTION_DOMAIN_CTX(ctx)); }


void
_gr_fraction_init(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    gr_init(a, domain_ctx);
    gr_init(b, domain_ctx);
    GR_MUST_SUCCEED(gr_one(b, domain_ctx));
}

void
_gr_fraction_clear(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    gr_clear(a, domain_ctx);
    gr_clear(b, domain_ctx);
}

int
_gr_fraction_write(gr_stream_t out, gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_stream_write(out, "(");
    status |= gr_write(out, a, domain_ctx);
    status |= gr_stream_write(out, ") / (");
    status |= gr_write(out, b, domain_ctx);
    status |= gr_stream_write(out, ")");

    return GR_SUCCESS;
}

void
_gr_fraction_set_shallow(gr_ptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_srcptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    gr_set_shallow(a, c, domain_ctx);
    gr_set_shallow(b, d, domain_ctx);
}

void
_gr_fraction_swap(gr_ptr x, gr_ptr y, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    gr_swap(a, c, domain_ctx);
    gr_swap(b, d, domain_ctx);
}

static int
gr_fraction_canonicalise(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (!(GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION))
    {
        gr_ptr t;
        GR_TMP_INIT(t, domain_ctx);

        status |= gr_gcd(t, a, b, domain_ctx);

        if (gr_is_one(t, domain_ctx) != T_TRUE)
        {
            status |= gr_divexact(a, a, t, domain_ctx);
            status |= gr_divexact(b, b, t, domain_ctx);
        }

        status |= gr_canonical_associate(b, t, b, domain_ctx);

        if (gr_is_one(t, domain_ctx) != T_TRUE)
            status |= gr_mul(a, a, t, domain_ctx);

        GR_TMP_CLEAR(t, domain_ctx);
    }

    return status;
}

static int
gr_fraction_canonicalise_unit(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (!(GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION))
    {
        gr_ptr t;
        GR_TMP_INIT(t, domain_ctx);

        status |= gr_canonical_associate(b, t, b, domain_ctx);

        if (gr_is_one(t, domain_ctx) != T_TRUE)
            status |= gr_mul(a, a, t, domain_ctx);

        GR_TMP_CLEAR(t, domain_ctx);
    }

    return status;
}


int
_gr_fraction_randtest(gr_ptr x, flint_rand_t state, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_randtest(a, state, domain_ctx);

    if (n_randint(state, 2))
    {
        status |= gr_one(b, domain_ctx);
    }
    else
    {
        status |= gr_randtest_not_zero(b, state, domain_ctx);
        status |= gr_fraction_canonicalise(x, ctx);
    }

    return status;
}

int
_gr_fraction_zero(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_zero(a, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_one(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_one(a, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_neg_one(gr_ptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg_one(a, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_set(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_srcptr c = NUMER(x, ctx), d = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(a, c, domain_ctx);
    status |= gr_set(b, d, domain_ctx);
    return status;
}

int
_gr_fraction_set_si(gr_ptr res, slong x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_si(a, x, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_set_ui(gr_ptr res, ulong x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_ui(a, x, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_set_fmpz(gr_ptr res, const fmpz_t x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_fmpz(a, x, domain_ctx);
    status |= gr_one(b, domain_ctx);
    return status;
}

int
_gr_fraction_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (x_ctx == domain_ctx)
    {
        status |= gr_set(a, x, domain_ctx);
        status |= gr_one(b, domain_ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_FRACTION)
    {
        gr_srcptr c = NUMER(x, x_ctx), d = DENOM(x, x_ctx);
        gr_ctx_struct * x_domain_ctx = DOMAIN(x_ctx);

        status |= gr_set_other(a, c, x_domain_ctx, domain_ctx);
        status |= gr_set_other(b, d, x_domain_ctx, domain_ctx);

        if (status == GR_SUCCESS)
            status = gr_fraction_canonicalise(res, ctx);
    }
    else
    {
        status = gr_set_other(a, x, x_ctx, domain_ctx);

        if (status == GR_SUCCESS)
        {
            status = gr_one(b, domain_ctx);
        }
        else
        {
            /* To do: recognize various other types with internal fractions
               more directly, e.g. nf, fmpq_poly, fmpq_mpoly, Frac(R)[x], ...?
            gr_ptr t, u;
            GR_TMP_INIT2(t, u, x_ctx);

            status = gr_numerator(t, x, x_ctx);
            status |= gr_denominator(u, x, x_ctx);

            if (status == GR_SUCCESS && gr_one(u, x_ctx) != T_TRUE)
            {
                status = gr_set_other(a, t, x_ctx, ctx);
                status |= gr_set_other(b, u, x_ctx, ctx);
                status |= gr_fraction_canonicalise(res, ctx);
            }
            else
                status = GR_UNABLE;

            GR_TMP_CLEAR2(t, u, x_ctx);

            if (status != GR_SUCCESS)
                status = gr_generic_set_other(res, x, x_ctx, ctx);
             */

            status = gr_generic_set_other(res, x, x_ctx, ctx);
        }
    }

    return status;
}

int
_gr_fraction_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    /* Get generators of the domain ring */
    gr_vec_init(vec1, 0, DOMAIN(ctx));
    status = gr_gens_recursive(vec1, DOMAIN(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n, ctx);


    /* Promote to fractions */
    for (i = 0; i < n; i++)
    {
        gr_ptr x = gr_vec_entry_ptr(vec, i, ctx);
        gr_srcptr y = gr_vec_entry_srcptr(vec1, i, DOMAIN(ctx));

        status |= gr_set(NUMER(x, ctx), y, DOMAIN(ctx));
        status |= gr_one(DENOM(x, ctx), DOMAIN(ctx));
    }

    gr_vec_clear(vec1, POLYNOMIAL_ELEM_CTX(ctx));

    return status;
}

truth_t
_gr_fraction_equal(gr_srcptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a, b, c, d;
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    a = NUMER(x, ctx);
    b = DENOM(x, ctx);
    c = NUMER(y, ctx);
    d = DENOM(y, ctx);

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_STRONGLY_CANONICAL)
    {
        truth_t eq1, eq2;

        eq1 = gr_equal(a, c, domain_ctx);
        if (eq1 == T_FALSE)
            return eq1;

        eq2 = gr_equal(b, d, domain_ctx);

        return truth_and(eq1, eq2);
    }
    else
    {
        /* ad = bc */
        gr_ptr t, u;
        truth_t eq;
        int status;

        GR_TMP_INIT2(t, u, domain_ctx);

        status = gr_mul(t, a, d, domain_ctx);
        status |= gr_mul(u, b, c, domain_ctx);

        eq = (status == GR_UNABLE) ? T_UNKNOWN : gr_equal(t, u, domain_ctx);

        GR_TMP_CLEAR2(t, u, domain_ctx);

        return eq;
    }
}

truth_t
_gr_fraction_is_zero(gr_srcptr x, gr_fraction_ctx_t ctx)
{
    return gr_is_zero(NUMER(x, ctx), DOMAIN(ctx));
}

truth_t
_gr_fraction_is_one(gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    /* XXX: this branch is not necessarily faster */
    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_STRONGLY_CANONICAL)
    {
        truth_t eq1, eq2;

        eq1 = gr_is_one(a, domain_ctx);
        if (eq1 == T_FALSE)
            return eq1;

        eq2 = gr_is_one(b, domain_ctx);

        return truth_and(eq1, eq2);
    }
    else
    {
        return gr_equal(a, b, domain_ctx);
    }
}

truth_t
_gr_fraction_is_neg_one(gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_STRONGLY_CANONICAL)
    {
        truth_t eq1, eq2;

        eq1 = gr_is_neg_one(a, domain_ctx);
        if (eq1 == T_FALSE)
            return eq1;

        eq2 = gr_is_one(b, domain_ctx);

        return truth_and(eq1, eq2);
    }
    else
    {
        /* a = -b */
        gr_ptr t;
        truth_t eq;
        int status;

        GR_TMP_INIT(t, domain_ctx);
        status = gr_neg(t, b, domain_ctx);
        eq = (status == GR_UNABLE) ? T_UNKNOWN : gr_equal(a, t, domain_ctx);
        GR_TMP_CLEAR(t, domain_ctx);

        return eq;
    }
}


int
_gr_fraction_neg(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_ptr a = NUMER(res, ctx), b = DENOM(res, ctx);
    gr_srcptr c = NUMER(x, ctx), d = DENOM(x, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg(a, c, domain_ctx);
    status |= gr_set(b, d, domain_ctx);
    return status;
}

int
_gr_fraction_add_early_reduction(gr_ptr res_num, gr_ptr res_den,
                gr_srcptr x_num, gr_srcptr x_den,
                gr_srcptr y_num, gr_srcptr y_den, int subtract, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (gr_is_zero(x_num, ctx) == T_TRUE)
    {
        status |= subtract ? gr_neg(res_num, y_num, ctx)
                           : gr_set(res_num, y_num, ctx);
        status |= gr_set(res_den, y_den, ctx);
        return status;
    }

    if (gr_is_zero(y_num, ctx) == T_TRUE)
    {
        status |= gr_set(res_num, x_num, ctx);
        status |= gr_set(res_den, x_den, ctx);
        return status;
    }

    /* XXX: is this special case worth it? */
    if (gr_equal(x_den, y_den, ctx) == T_TRUE)
    {
        status |= subtract ? gr_sub(res_num, x_num, y_num, ctx)
                           : gr_add(res_num, x_num, y_num, ctx);

        if (gr_is_one(x_den, ctx) == T_TRUE || gr_is_zero(res_num, ctx) == T_TRUE)
        {
            status |= gr_one(res_den, ctx);
        }
        else
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);

            status |= gr_gcd(t, res_num, x_den, ctx);

            if (gr_is_one(t, ctx) == T_TRUE)
            {
                status |= gr_set(res_den, x_den, ctx);
            }
            else
            {
                status |= gr_divexact(res_num, res_num, t, ctx);
                status |= gr_divexact(res_den, x_den, t, ctx);
            }

            status |= gr_canonical_associate(res_den, t, res_den, ctx);
            if (gr_is_one(t, ctx) != T_TRUE)
                status |= gr_mul(res_num, res_num, t, ctx);

            GR_TMP_CLEAR(t, ctx);
        }

        return status;
    }

    if (gr_is_one(x_den, ctx) == T_TRUE)
    {
        if (res_num == y_num)
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);

            status |= gr_mul(t, x_num, y_den, ctx);
            status |= subtract ? gr_sub(res_num, t, y_num, ctx)
                               : gr_add(res_num, t, y_num, ctx);

            GR_TMP_CLEAR(t, ctx);
        }
        else
        {
            status |= gr_mul(res_num, x_num, y_den, ctx);
            status |= subtract ? gr_sub(res_num, res_num, y_num, ctx)
                               : gr_add(res_num, res_num, y_num, ctx);
        }

        if (res_den != y_den)
            status |= gr_set(res_den, y_den, ctx);
        return status;
    }

    if (gr_is_one(y_den, ctx) == T_TRUE)
    {
        if (res_num == x_num)
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);
            status |= gr_mul(t, x_den, y_num, ctx);
            status |= subtract ? gr_sub(res_num, x_num, t, ctx)
                               : gr_add(res_num, x_num, t, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else
        {
            status |= gr_mul(res_num, x_den, y_num, ctx);
            status |= subtract ? gr_sub(res_num, x_num, res_num, ctx)
                               : gr_add(res_num, x_num, res_num, ctx);
        }

        if (res_den != x_den)
            status |= gr_set(res_den, x_den, ctx);
        return status;
    }

    {
        gr_ptr g;
        GR_TMP_INIT(g, ctx);

        status |= gr_gcd(g, x_den, y_den, ctx);

        if (gr_is_one(g, ctx) == T_TRUE)
        {
            gr_ptr t, u;
            GR_TMP_INIT2(t, u, ctx);

            /* todo: avoid one alloc */
            status |= gr_mul(t, x_num, y_den, ctx);
            status |= gr_mul(u, y_num, x_den, ctx);
            status |= subtract ? gr_sub(res_num, t, u, ctx)
                               : gr_add(res_num, t, u, ctx);
            status |= gr_mul(res_den, x_den, y_den, ctx);

            GR_TMP_CLEAR2(t, u, ctx);
        }
        else
        {
            gr_ptr a, b, t, u;
            GR_TMP_INIT4(a, b, t, u, ctx);

            status |= gr_divexact(a, x_den, g, ctx);
            status |= gr_divexact(b, y_den, g, ctx);

            status |= gr_mul(t, x_num, b, ctx);
            status |= gr_mul(u, y_num, a, ctx);
            status |= subtract ? gr_sub(res_num, t, u, ctx)
                               : gr_add(res_num, t, u, ctx);

            status |= gr_gcd(t, res_num, g, ctx);

            if (gr_is_one(t, ctx) == T_TRUE)
            {
                status |= gr_mul(res_den, x_den, b, ctx);
            }
            else
            {
                status |= gr_divexact(res_num, res_num, t, ctx);
                status |= gr_divexact(g, x_den, t, ctx);
                status |= gr_mul(res_den, g, b, ctx);
            }

            GR_TMP_CLEAR4(a, b, t, u, ctx);
        }

        status |= gr_canonical_associate(res_den, g, res_den, ctx);
        if (gr_is_one(g, ctx) != T_TRUE)
            status |= gr_mul(res_num, res_num, g, ctx);

        GR_TMP_CLEAR(g, ctx);

        return status;
    }
}

int
_gr_fraction_mul_early_reduction(gr_ptr res_num, gr_ptr res_den,
                gr_srcptr x_num, gr_srcptr x_den,
                gr_srcptr y_num, gr_srcptr y_den, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (gr_is_zero(x_num, ctx) == T_TRUE || gr_is_zero(y_num, ctx) == T_TRUE)
    {
        status |= gr_zero(res_num, ctx);
        status |= gr_one(res_den, ctx);
        return status;
    }

    /* Also covers squaring */
    if (x_den == y_den || gr_equal(x_den, y_den, ctx) == T_TRUE)
    {
        status |= gr_mul(res_num, x_num, y_num, ctx);
        status |= gr_mul(res_den, x_den, y_den, ctx);

        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status |= gr_canonical_associate(res_den, t, res_den, ctx);
        if (gr_is_one(t, ctx) != T_TRUE)
            status |= gr_mul(res_num, res_num, t, ctx);
        GR_TMP_CLEAR(t, ctx);

        return status;
    }

    if (gr_is_one(x_den, ctx) == T_TRUE)
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, ctx);

        status |= gr_gcd(t, x_num, y_den, ctx);

        if (gr_is_one(t, ctx) == T_TRUE)
        {
            status |= gr_mul(res_num, x_num, y_num, ctx);
            status |= gr_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            status |= gr_divexact(u, x_num, t, ctx);
            status |= gr_mul(res_num, u, y_num, ctx);
            status |= gr_divexact(u, y_den, t, ctx);
            status |= gr_mul(res_den, x_den, u, ctx);
        }

        status |= gr_canonical_associate(res_den, t, res_den, ctx);
        if (gr_is_one(t, ctx) != T_TRUE)
            status |= gr_mul(res_num, res_num, t, ctx);

        GR_TMP_CLEAR2(t, u, ctx);
        return status;
    }

    if (gr_is_one(y_den, ctx) == T_TRUE)
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, ctx);

        status |= gr_gcd(t, y_num, x_den, ctx);

        if (gr_is_one(t, ctx) == T_TRUE)
        {
            status |= gr_mul(res_num, x_num, y_num, ctx);
            status |= gr_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            status |= gr_divexact(u, y_num, t, ctx);
            status |= gr_mul(res_num, u, x_num, ctx);
            status |= gr_divexact(u, x_den, t, ctx);
            status |= gr_mul(res_den, y_den, u, ctx);
        }

        status |= gr_canonical_associate(res_den, t, res_den, ctx);
        if (gr_is_one(t, ctx) != T_TRUE)
            status |= gr_mul(res_num, res_num, t, ctx);

        GR_TMP_CLEAR2(t, u, ctx);
        return status;
    }


    {
        gr_ptr t, u, x, y;

        GR_TMP_INIT4(t, u, x, y, ctx);

        status |= gr_gcd(t, x_num, y_den, ctx);

        if (gr_is_one(t, ctx) == T_TRUE)
        {
            status |= gr_gcd(u, x_den, y_num, ctx);

            if (gr_is_one(u, ctx) == T_TRUE)
            {
                status |= gr_mul(res_num, x_num, y_num, ctx);
                status |= gr_mul(res_den, x_den, y_den, ctx);
            }
            else
            {
                status |= gr_divexact(y, y_num, u, ctx);
                status |= gr_mul(res_num, x_num, y, ctx);
                status |= gr_divexact(x, x_den, u, ctx);
                status |= gr_mul(res_den, x, y_den, ctx);
            }
        }
        else
        {
            status |= gr_gcd(u, x_den, y_num, ctx);

            if (gr_is_one(u, ctx) == T_TRUE)
            {
                status |= gr_divexact(x, x_num, t, ctx);
                status |= gr_mul(res_num, x, y_num, ctx);
                status |= gr_divexact(y, y_den, t, ctx);
                status |= gr_mul(res_den, x_den, y, ctx);
            }
            else
            {
                status |= gr_divexact(x, x_num, t, ctx);
                status |= gr_divexact(y, y_num, u, ctx);
                status |= gr_mul(res_num, x, y, ctx);
                status |= gr_divexact(x, x_den, u, ctx);
                status |= gr_divexact(y, y_den, t, ctx);
                status |= gr_mul(res_den, x, y, ctx);
            }
        }

        status |= gr_canonical_associate(res_den, t, res_den, ctx);
        if (gr_is_one(t, ctx) != T_TRUE)
            status |= gr_mul(res_num, res_num, t, ctx);

        GR_TMP_CLEAR4(t, u, x, y, ctx);
        return status;
    }
}

int
_gr_fraction_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_srcptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION)
    {
        /* p/q = (ad + bc) / bd */
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, domain_ctx);
        status |= gr_mul(t, a, d, domain_ctx);
        status |= gr_mul(u, b, c, domain_ctx);
        status |= gr_mul(q, b, d, domain_ctx);
        status |= gr_add(p, t, u, domain_ctx);
        GR_TMP_CLEAR2(t, u, domain_ctx);
    }
    else
    {
        status = _gr_fraction_add_early_reduction(p, q,
            a, b, c, d, 0, domain_ctx);
    }

    return status;
}

int
_gr_fraction_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_srcptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION)
    {
        /* p/q = (ad - bc) / bd */
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, domain_ctx);
        status |= gr_mul(t, a, d, domain_ctx);
        status |= gr_mul(u, b, c, domain_ctx);
        status |= gr_mul(q, b, d, domain_ctx);
        status |= gr_sub(p, t, u, domain_ctx);
        GR_TMP_CLEAR2(t, u, domain_ctx);
    }
    else
    {
        status = _gr_fraction_add_early_reduction(p, q,
            a, b, c, d, 1, domain_ctx);
    }

    return status;
}

int
_gr_fraction_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_srcptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION)
    {
        status |= gr_mul(p, a, c, domain_ctx);
        status |= gr_mul(q, b, d, domain_ctx);
    }
    else
    {
        status = _gr_fraction_mul_early_reduction(p, q,
            a, b, c, d, domain_ctx);
    }

    return status;
}

int
_gr_fraction_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_srcptr c = NUMER(y, ctx), d = DENOM(y, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    truth_t zero;

    zero = gr_is_zero(c, domain_ctx);

    if (zero == T_TRUE)
        return GR_DOMAIN;
    if (zero == T_UNKNOWN)
        return GR_UNABLE;

    if (GR_FRACTION_FLAGS(ctx) & GR_FRACTION_NO_REDUCTION)
    {
        if (res == x || res == y)
        {
            gr_ptr t;
            GR_TMP_INIT(t, domain_ctx);
            status |= gr_mul(t, a, d, domain_ctx);
            status |= gr_mul(q, b, c, domain_ctx);
            gr_swap(p, t, domain_ctx);
            GR_TMP_CLEAR(t, domain_ctx);
        }
        else
        {
            status |= gr_mul(p, a, d, domain_ctx);
            status |= gr_mul(q, b, c, domain_ctx);
        }
    }
    else
    {
        if (res == y)
        {
            gr_ptr t, u;
            GR_TMP_INIT2(t, u, domain_ctx);
            status |= gr_set(t, c, domain_ctx);
            status |= gr_set(u, d, domain_ctx);
            status |= _gr_fraction_mul_early_reduction(p, q,
                    a, b, u, t, domain_ctx);
            GR_TMP_CLEAR2(t, u, domain_ctx);
        }
        else
        {
            status = _gr_fraction_mul_early_reduction(p, q,
                a, b, d, c, domain_ctx);
        }
    }

    return status;
}

int
_gr_fraction_inv(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    truth_t zero;

    zero = gr_is_zero(a, domain_ctx);

    if (zero == T_TRUE)
        return GR_DOMAIN;
    if (zero == T_UNKNOWN)
        return GR_UNABLE;

    if (res == x)
    {
        gr_swap(p, q, domain_ctx);
    }
    else
    {
        status |= gr_set(p, b, domain_ctx);
        status |= gr_set(q, a, domain_ctx);
    }

    status |= gr_fraction_canonicalise_unit(res, ctx);

    return status;
}

truth_t
_gr_fraction_is_invertible(gr_srcptr x, gr_fraction_ctx_t ctx)
{
    return truth_not(gr_is_zero(NUMER(x, ctx), DOMAIN(ctx)));
}

int
_gr_fraction_sqrt(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_sqrt(p, a, domain_ctx);
    status |= gr_sqrt(q, b, domain_ctx);
    status |= gr_fraction_canonicalise_unit(res, ctx);

    /* todo: when can we guarantee GR_DOMAIN? */
    if (status != GR_SUCCESS)
    {
        /* don't accidentally construct the invalid fraction 1 / 0 */
        GR_IGNORE(gr_zero(res, ctx));
        status = GR_UNABLE;
    }

    return status;
}

int
_gr_fraction_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    status |= gr_pow_ui(p, a, y, domain_ctx);
    status |= gr_pow_ui(q, b, y, domain_ctx);
    status |= gr_fraction_canonicalise_unit(res, ctx);
    return status;
}

int
_gr_fraction_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_fraction_ctx_t ctx)
{
    gr_srcptr a = NUMER(x, ctx), b = DENOM(x, ctx);
    gr_ptr p = NUMER(res, ctx), q = DENOM(res, ctx);
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);
    int status = GR_SUCCESS;

    if (fmpz_sgn(y) >= 0)
    {
        status |= gr_pow_fmpz(p, a, y, domain_ctx);
        status |= gr_pow_fmpz(q, b, y, domain_ctx);
    }
    else
    {
        truth_t zero;
        fmpz_t e;

        zero = gr_is_zero(a, domain_ctx);
        if (zero == T_TRUE)
            return GR_DOMAIN;
        if (zero == T_UNKNOWN)
            return GR_UNABLE;

        fmpz_init(e);
        fmpz_neg(e, y);

        status |= gr_pow_fmpz(p, a, e, domain_ctx);
        status |= gr_pow_fmpz(q, b, e, domain_ctx);
        gr_swap(p, q, domain_ctx);
        fmpz_clear(e);
    }

    status |= gr_fraction_canonicalise_unit(res, ctx);
    return status;
}

int
_gr_fraction_numerator(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    status |= gr_set(NUMER(res, ctx), NUMER(x, ctx), domain_ctx);
    status |= gr_one(DENOM(res, ctx), domain_ctx);
    return status;
}

int
_gr_fraction_denominator(gr_ptr res, gr_srcptr x, gr_fraction_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    status |= gr_set(NUMER(res, ctx), DENOM(x, ctx), domain_ctx);
    status |= gr_one(DENOM(res, ctx), domain_ctx);
    return status;
}

int
_gr_fraction_i(gr_ptr res, gr_fraction_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    status |= gr_i(NUMER(res, ctx), domain_ctx);
    status |= gr_one(DENOM(res, ctx), domain_ctx);

    return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_fraction_pi(gr_ptr res, gr_fraction_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * domain_ctx = DOMAIN(ctx);

    status |= gr_pi(NUMER(res, ctx), domain_ctx);
    status |= gr_one(DENOM(res, ctx), domain_ctx);

    return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
}

int _gr_fraction_methods_initialized = 0;

gr_static_method_table _gr_fraction_methods;

gr_method_tab_input _gr_fraction_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fraction_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fraction_ctx_clear},


    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_fraction_ctx_is_certainly_field_else_unknown},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr)  _gr_fraction_ctx_is_certainly_field_else_unknown},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr)  _gr_fraction_ctx_is_certainly_field_else_unknown},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,  (gr_funcptr)  _gr_fraction_ctx_is_certainly_field_else_unknown},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr)  _gr_fraction_ctx_is_certainly_field_else_unknown},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr)  _gr_fraction_ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) _gr_fraction_ctx_is_finite},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) _gr_fraction_ctx_is_finite_characteristic},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) _gr_fraction_ctx_is_exact},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_fraction_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fraction_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fraction_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fraction_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fraction_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fraction_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fraction_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fraction_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_fraction_neg_one},
/*
    todo: should this inherit the base ring's gen()?
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fraction_gen},
*/
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) _gr_fraction_gens_recursive},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fraction_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fraction_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fraction_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fraction_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fraction_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fraction_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fraction_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fraction_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fraction_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
/*
    {GR_METHOD_SET_FEXPR,       (gr_funcptr) _gr_fraction_set_fexpr},
    {GR_METHOD_GET_FEXPR,       (gr_funcptr) _gr_fraction_get_fexpr},
*/
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fraction_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_fraction_add},
/*
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_fraction_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_fraction_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_fraction_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_fraction_add_fmpq},
*/
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fraction_sub},
/*
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_fraction_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_fraction_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_fraction_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_fraction_sub_fmpq},
*/
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fraction_mul},
/*
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fraction_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fraction_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fraction_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_fraction_mul_fmpq},

    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fraction_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fraction_sqr},
*/

    {GR_METHOD_DIV,             (gr_funcptr) _gr_fraction_div},
/*
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_fraction_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_fraction_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_fraction_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_fraction_div_fmpq},
*/

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fraction_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fraction_inv},

    {GR_METHOD_SQRT,             (gr_funcptr) _gr_fraction_sqrt},

    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fraction_pow_ui},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fraction_pow_fmpz},

    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_fraction_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_fraction_denominator},

    {GR_METHOD_I,               (gr_funcptr) _gr_fraction_i},
    {GR_METHOD_PI,              (gr_funcptr) _gr_fraction_pi},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_fraction(gr_ctx_t ctx, gr_ctx_t domain, int flags)
{
    ctx->which_ring = GR_CTX_GR_FRACTION;
    ctx->sizeof_elem = 2 * domain->sizeof_elem;
    ctx->size_limit = WORD_MAX;

    GR_FRACTION_DOMAIN_CTX(ctx) = domain;
    GR_FRACTION_FLAGS(ctx) = flags;

    ctx->methods = _gr_fraction_methods;

    if (!_gr_fraction_methods_initialized)
    {
        gr_method_tab_init(_gr_fraction_methods, _gr_fraction_methods_input);
        _gr_fraction_methods_initialized = 1;
    }
}

