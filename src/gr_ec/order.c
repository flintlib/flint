/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The order of E(F_q), remembered in the context.

    Point counting is expensive and its answer is used by almost everything
    else -- reducing a scalar before a multiplication, deciding whether a
    division by n is unique, accepting a ring of scalars -- so the context
    caches it. Two strengths are distinguished:

      GR_EC_ORDER_EXACT         the value is #E(F_q)
      GR_EC_ORDER_ANNIHILATOR   the value is some multiple of the group
                                exponent

    The weaker one is what arithmetic actually needs: if m kills every point
    then k P = (k mod m) P, and n is invertible on the group as soon as
    gcd(n, m) = 1. It is also the most that can be established without
    counting, which matters because a caller who hands over a modulus should
    not thereby trigger a point count.

    Nothing here is reached by the counting functions themselves:
    gr_ec_ctx_cardinality and friends always do the work. gr_ec_ctx_order is
    the cached entry point, and is what the generic interface calls.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/*
    How many random points a claimed annihilator is tested against.

    The points killed by m form a subgroup, so if m does not kill the whole
    group it misses at least half of it and a random point refutes m with
    probability at least 1/2. Twenty of them leave a chance below 10^-6,
    and cost twenty scalar multiplications -- against a point count, which
    is what the alternative would be, that is cheap at any size worth
    caring about.
*/
#define GR_EC_ORDER_PROBE_POINTS 20

/*
    Below this the field can simply be walked, so a claimed order can be
    checked against the truth instead of against random points. Matches the
    cutoff gr_ec_ctx_cardinality uses for the same reason.
*/
#define GR_EC_ORDER_EXACT_CHECK_MAX_Q WORD(10000)

void
gr_ec_ctx_clear_order(gr_ec_ctx_t ctx)
{
    fmpz_zero(&GR_EC_CTX(ctx)->order);
    GR_EC_CTX(ctx)->order_kind = GR_EC_ORDER_UNKNOWN;
}

gr_ec_order_kind_t
gr_ec_ctx_get_cached_order(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ec_order_kind_t kind = GR_EC_CTX(ctx)->order_kind;

    if (kind != GR_EC_ORDER_UNKNOWN)
        fmpz_set(res, &GR_EC_CTX(ctx)->order);

    return kind;
}

/*
    Record what has been learned, without ever weakening what is already
    there. Two annihilators combine: a Bezout combination of them kills
    every point that both do, so their gcd is an annihilator as well, and a
    smaller one is a better scalar reduction.
*/
static void
_gr_ec_ctx_store(gr_ec_ctx_t ctx, const fmpz_t value, gr_ec_order_kind_t kind)
{
    gr_ec_order_kind_t have = GR_EC_CTX(ctx)->order_kind;

    if (kind == GR_EC_ORDER_EXACT)
    {
        fmpz_set(&GR_EC_CTX(ctx)->order, value);
        GR_EC_CTX(ctx)->order_kind = GR_EC_ORDER_EXACT;
        return;
    }

    if (have == GR_EC_ORDER_EXACT)
        return;

    if (have == GR_EC_ORDER_ANNIHILATOR)
        fmpz_gcd(&GR_EC_CTX(ctx)->order, &GR_EC_CTX(ctx)->order, value);
    else
        fmpz_set(&GR_EC_CTX(ctx)->order, value);

    GR_EC_CTX(ctx)->order_kind = GR_EC_ORDER_ANNIHILATOR;
}

/*
    Does m kill GR_EC_ORDER_PROBE_POINTS random points?

    Returns 1 if it kills all of them, 0 if some point survives, and -1 if
    no point could be produced at all -- which happens over a base ring
    without square roots, since a random point needs one.
*/
static int
_gr_ec_probe_annihilator(const fmpz_t m, gr_ec_ctx_t ctx)
{
    flint_rand_t state;
    gr_ec_point_t P, T;
    slong tries, tested = 0;
    int result = -1;

    flint_rand_init(state);
    flint_rand_set_seed(state, UWORD(0x6a09e667f3bcc908),
            UWORD(0xbb67ae8584caa73b));

    gr_ec_point_init(P, ctx);
    gr_ec_point_init(T, ctx);

    for (tries = 0; tries < 4 * GR_EC_ORDER_PROBE_POINTS
            && tested < GR_EC_ORDER_PROBE_POINTS; tries++)
    {
        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS)
            break;

        /* the point at infinity says nothing */
        if (gr_ec_point_is_inf(P, ctx) != T_FALSE)
            continue;

        if (gr_ec_point_mul_fmpz(T, P, m, ctx) != GR_SUCCESS)
            break;

        tested++;

        if (gr_ec_point_is_inf(T, ctx) != T_TRUE)
        {
            result = 0;
            goto cleanup;
        }
    }

    if (tested > 0)
        result = 1;

cleanup:
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(P, ctx);
    flint_rand_clear(state);

    return result;
}

truth_t
gr_ec_ctx_annihilates(gr_ec_ctx_t ctx, const fmpz_t m)
{
    int probe;

    if (fmpz_sgn(m) <= 0)
        return T_FALSE;

    /* a multiple of something that kills the group kills it too */
    if (GR_EC_CTX(ctx)->order_kind != GR_EC_ORDER_UNKNOWN
            && fmpz_divisible(m, &GR_EC_CTX(ctx)->order))
        return T_TRUE;

    probe = _gr_ec_probe_annihilator(m, ctx);

    if (probe < 0)
        return T_UNKNOWN;

    if (probe == 0)
        return T_FALSE;

    /* evidence about m, not a point count, so it goes in as the weaker
       kind and never displaces a count */
    _gr_ec_ctx_store(ctx, m, GR_EC_ORDER_ANNIHILATOR);

    return T_TRUE;
}

int
gr_ec_ctx_set_annihilator(gr_ec_ctx_t ctx, const fmpz_t m)
{
    truth_t ok = gr_ec_ctx_annihilates(ctx, m);

    if (ok == T_FALSE)
        return GR_DOMAIN;

    /* unverifiable is still usable: the caller asserted it */
    _gr_ec_ctx_store(ctx, m, GR_EC_ORDER_ANNIHILATOR);

    return GR_SUCCESS;
}

/* |N - (q + 1)| <= 2 sqrt(q) */
static int
_gr_ec_hasse_ok(const fmpz_t N, const fmpz_t q)
{
    fmpz_t t, bound;
    int ok;

    fmpz_init(t);
    fmpz_init(bound);

    fmpz_sub(t, N, q);
    fmpz_sub_ui(t, t, 1);
    fmpz_mul(t, t, t);
    fmpz_mul_ui(bound, q, 4);

    ok = (fmpz_cmp(t, bound) <= 0);

    fmpz_clear(t);
    fmpz_clear(bound);

    return ok;
}

int
gr_ec_ctx_set_order(gr_ec_ctx_t ctx, const fmpz_t N)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q;
    int status = GR_SUCCESS;

    if (fmpz_sgn(N) <= 0)
        return GR_DOMAIN;

    /* already counted: the claim is either right or wrong, no probing */
    if (GR_EC_CTX(ctx)->order_kind == GR_EC_ORDER_EXACT)
        return fmpz_equal(N, &GR_EC_CTX(ctx)->order) ? GR_SUCCESS : GR_DOMAIN;

    fmpz_init(q);

    if (gr_ctx_cardinality_fmpz(q, R) == GR_SUCCESS)
    {
        if (!_gr_ec_hasse_ok(N, q))
        {
            status = GR_DOMAIN;
            goto cleanup;
        }

        /* small enough to settle outright */
        if (fmpz_cmp_si(q, GR_EC_ORDER_EXACT_CHECK_MAX_Q) <= 0)
        {
            fmpz_t truth;
            fmpz_init(truth);

            if (gr_ec_ctx_cardinality_naive(truth, ctx) == GR_SUCCESS)
            {
                status = fmpz_equal(N, truth) ? GR_SUCCESS : GR_DOMAIN;

                if (status == GR_SUCCESS)
                    _gr_ec_ctx_store(ctx, truth, GR_EC_ORDER_EXACT);

                fmpz_clear(truth);
                goto cleanup;
            }

            fmpz_clear(truth);
        }
    }

    /* otherwise check as far as is cheap: N must kill the group */
    if (_gr_ec_probe_annihilator(N, ctx) == 0)
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    _gr_ec_ctx_store(ctx, N, GR_EC_ORDER_EXACT);

cleanup:
    fmpz_clear(q);

    return status;
}

int
gr_ec_ctx_order(fmpz_t res, gr_ec_ctx_t ctx)
{
    int status;

    if (GR_EC_CTX(ctx)->order_kind == GR_EC_ORDER_EXACT)
    {
        fmpz_set(res, &GR_EC_CTX(ctx)->order);
        return GR_SUCCESS;
    }

    status = gr_ec_ctx_cardinality(res, ctx);

    if (status == GR_SUCCESS)
    {
        /* a count beats anything a caller asserted, so store it outright */
        fmpz_set(&GR_EC_CTX(ctx)->order, res);
        GR_EC_CTX(ctx)->order_kind = GR_EC_ORDER_EXACT;
    }

    return status;
}
