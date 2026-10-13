/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* number of random points a claimed cardinality is tested against */
#define GR_EC_CARDINALITY_PROBE_POINTS 20

/* below this, a claimed cardinality is checked by counting naively */
#define GR_EC_CARDINALITY_CHECK_MAX_Q WORD(10000)

/*
    Returns 0 if some random point is not killed by m, 1 if all the points
    tested are, and -1 if no point could be produced.
*/
static int
_gr_ec_probe_annihilator(const fmpz_t m, gr_ec_ctx_t ctx)
{
    flint_rand_t state;
    gr_ec_jac_point_t P, T;
    slong tries, tested = 0;
    int result = -1;

    flint_rand_init(state);

    gr_ec_jac_point_init(P, ctx);
    gr_ec_jac_point_init(T, ctx);

    for (tries = 0; tries < 4 * GR_EC_CARDINALITY_PROBE_POINTS
            && tested < GR_EC_CARDINALITY_PROBE_POINTS; tries++)
    {
        if (gr_ec_jac_point_randtest(P, state, ctx) != GR_SUCCESS)
            break;

        if (gr_ec_jac_point_is_inf(P, ctx) != T_FALSE)
            continue;

        if (gr_ec_jac_point_mul_fmpz(T, P, m, ctx) != GR_SUCCESS)
            break;

        tested++;

        if (gr_ec_jac_point_is_inf(T, ctx) != T_TRUE)
        {
            result = 0;
            break;
        }
    }

    if (result != 0 && tested > 0)
        result = 1;

    gr_ec_jac_point_clear(T, ctx);
    gr_ec_jac_point_clear(P, ctx);
    flint_rand_clear(state);

    return result;
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
gr_ec_ctx_set_cardinality(gr_ec_ctx_t ctx, const fmpz_t N)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q, M;
    int status = GR_SUCCESS;

    if (fmpz_sgn(N) <= 0)
        return GR_DOMAIN;

    if (!fmpz_is_zero(&GR_EC_CTX(ctx)->cardinality))
        return fmpz_equal(N, &GR_EC_CTX(ctx)->cardinality) ? GR_SUCCESS : GR_DOMAIN;

    fmpz_init(q);
    fmpz_init(M);

    if (gr_ctx_is_field(R) == T_TRUE && gr_ctx_cardinality_fmpz(q, R) == GR_SUCCESS)
    {
        if (!_gr_ec_hasse_ok(N, q))
            status = GR_DOMAIN;
        else if (fmpz_cmp_si(q, GR_EC_CARDINALITY_CHECK_MAX_Q) <= 0
                && gr_ec_ctx_cardinality_naive(M, ctx) == GR_SUCCESS
                && !fmpz_equal(N, M))
            status = GR_DOMAIN;
    }

    if (status == GR_SUCCESS && _gr_ec_probe_annihilator(N, ctx) == 0)
        status = GR_DOMAIN;

    if (status == GR_SUCCESS)
        fmpz_set(&GR_EC_CTX(ctx)->cardinality, N);

    fmpz_clear(q);
    fmpz_clear(M);

    return status;
}
