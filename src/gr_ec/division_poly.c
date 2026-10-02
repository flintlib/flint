/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "impl.h"

/*
    Division polynomials.

    For the general Weierstrass equation psi_1 = 1, psi_2 = 2y + a1 x + a3,
    and psi_n is a polynomial in x alone when n is odd, and psi_2 times a
    polynomial in x when n is even. We always return the univariate

        Psi_n = psi_n            (n odd)
        Psi_n = psi_n / psi_2    (n even)

    so that the result lies in R[x] for every n. Writing

        W = psi_2^2 = 4 x^3 + b2 x^2 + 2 b4 x + b6

    the recursions become, with no division anywhere, so that they are
    valid over any commutative ring including characteristic 2 and 3:

        Psi_{2m}   = Psi_m (Psi_{m+2} Psi_{m-1}^2 - Psi_{m-2} Psi_{m+1}^2)
        Psi_{2m+1} = Psi_{m+2} Psi_m^3 - W^2 Psi_{m-1} Psi_{m+1}^3   (m odd)
                     W^2 Psi_{m+2} Psi_m^3 - Psi_{m-1} Psi_{m+1}^3   (m even)

    A single Psi_n needs only the transitive closure of these recursions,
    which is about ten indices per halving of n: 55 polynomials for
    n = 1000 rather than the 1001 of the full table, and 132 for n = 10^6.
    We compute the interval of indices touched at each halving top-down,
    mark exactly which of them are needed, then evaluate bottom-up and
    free each level once the level above it has been built. The cost is
    dominated by the last step, so it is O(M(n^2)) rather than the
    O(n M(n^2)) of building the table.
*/

/* enough for any n that could conceivably be computed */
#define DIVPOLY_MAX_LEVELS 64

/* state bits */
#define HAVE_PSI 1
#define HAVE_SQ 2
#define HAVE_CB 4
#define NEEDED 8

typedef struct
{
    gr_poly_struct * psi;
    gr_poly_struct * sq;
    gr_poly_struct * cb;
    unsigned char * state;
    slong off[DIVPOLY_MAX_LEVELS + 1];
    ulong lo[DIVPOLY_MAX_LEVELS];
    ulong hi[DIVPOLY_MAX_LEVELS];
    slong nlev;
    slong total;
}
divpoly_state_struct;

#define SLOT(S, j, k) ((S)->off[j] + (slong) ((k) - (S)->lo[j]))

/* res = 4 x^3 + b2 x^2 + 2 b4 x + b6 */
int
gr_ec_ctx_psi2_sqr(gr_poly_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr c;
    int status = GR_SUCCESS;

    gr_poly_fit_length(res, 4, R);
    _gr_poly_set_length(res, 4, R);
    c = res->coeffs;

    status |= gr_set(GR_ENTRY(c, 0, sz), GR_EC_B6(ctx), R);
    status |= gr_mul_two(GR_ENTRY(c, 1, sz), GR_EC_B4(ctx), R);
    status |= gr_set(GR_ENTRY(c, 2, sz), GR_EC_B2(ctx), R);
    status |= gr_set_ui(GR_ENTRY(c, 3, sz), 4, R);

    _gr_poly_normalise(res, R);

    return status;
}

/* Psi_n for n <= 4 */
static int
_divpoly_base(gr_poly_t res, ulong n, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_srcptr b2 = GR_EC_B2(ctx), b4 = GR_EC_B4(ctx);
    gr_srcptr b6 = GR_EC_B6(ctx), b8 = GR_EC_B8(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr c, t;
    int status = GR_SUCCESS;

    if (n == 0)
        return gr_poly_zero(res, R);

    if (n == 1 || n == 2)
        return gr_poly_one(res, R);

    if (n == 3)
    {
        /* 3 x^4 + b2 x^3 + 3 b4 x^2 + 3 b6 x + b8 */
        gr_poly_fit_length(res, 5, R);
        _gr_poly_set_length(res, 5, R);
        c = res->coeffs;

        status |= gr_set(GR_ENTRY(c, 0, sz), b8, R);
        status |= gr_mul_ui(GR_ENTRY(c, 1, sz), b6, 3, R);
        status |= gr_mul_ui(GR_ENTRY(c, 2, sz), b4, 3, R);
        status |= gr_set(GR_ENTRY(c, 3, sz), b2, R);
        status |= gr_set_ui(GR_ENTRY(c, 4, sz), 3, R);

        _gr_poly_normalise(res, R);
        return status;
    }

    /*
        2 x^6 + b2 x^5 + 5 b4 x^4 + 10 b6 x^3 + 10 b8 x^2
            + (b2 b8 - b4 b6) x + (b4 b8 - b6^2)
    */
    gr_poly_fit_length(res, 7, R);
    _gr_poly_set_length(res, 7, R);
    c = res->coeffs;

    GR_TMP_INIT(t, R);

    status |= gr_mul(GR_ENTRY(c, 0, sz), b4, b8, R);
    status |= gr_sqr(t, b6, R);
    status |= gr_sub(GR_ENTRY(c, 0, sz), GR_ENTRY(c, 0, sz), t, R);

    status |= gr_mul(GR_ENTRY(c, 1, sz), b2, b8, R);
    status |= gr_mul(t, b4, b6, R);
    status |= gr_sub(GR_ENTRY(c, 1, sz), GR_ENTRY(c, 1, sz), t, R);

    status |= gr_mul_ui(GR_ENTRY(c, 2, sz), b8, 10, R);
    status |= gr_mul_ui(GR_ENTRY(c, 3, sz), b6, 10, R);
    status |= gr_mul_ui(GR_ENTRY(c, 4, sz), b4, 5, R);
    status |= gr_set(GR_ENTRY(c, 5, sz), b2, R);
    status |= gr_set_ui(GR_ENTRY(c, 6, sz), 2, R);

    GR_TMP_CLEAR(t, R);

    _gr_poly_normalise(res, R);

    return status;
}

static gr_poly_struct *
_divpoly_sq(divpoly_state_struct * S, slong j, ulong k, gr_ctx_t R, int * status)
{
    slong s = SLOT(S, j, k);

    if (!(S->state[s] & HAVE_SQ))
    {
        *status |= gr_poly_mul(S->sq + s, S->psi + s, S->psi + s, R);
        S->state[s] |= HAVE_SQ;
    }

    return S->sq + s;
}

static gr_poly_struct *
_divpoly_cb(divpoly_state_struct * S, slong j, ulong k, gr_ctx_t R, int * status)
{
    slong s = SLOT(S, j, k);

    if (!(S->state[s] & HAVE_CB))
    {
        gr_poly_struct * q = _divpoly_sq(S, j, k, R, status);
        *status |= gr_poly_mul(S->cb + s, q, S->psi + s, R);
        S->state[s] |= HAVE_CB;
    }

    return S->cb + s;
}

/*
    Psi_k from the five (four for odd k) entries around k/2 at level j,
    which the caller has already built.
*/
static int
_divpoly_step(gr_poly_t res, divpoly_state_struct * S, slong j, ulong k,
        const gr_poly_t W2, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    ulong m = k / 2;
    gr_poly_t t, u;
    int status = GR_SUCCESS;

    gr_poly_init(t, R);
    gr_poly_init(u, R);

    if (k % 2 == 0)
    {
        /* Psi_m (Psi_{m+2} Psi_{m-1}^2 - Psi_{m-2} Psi_{m+1}^2) */
        status |= gr_poly_mul(t, S->psi + SLOT(S, j, m + 2),
                    _divpoly_sq(S, j, m - 1, R, &status), R);
        status |= gr_poly_mul(u, S->psi + SLOT(S, j, m - 2),
                    _divpoly_sq(S, j, m + 1, R, &status), R);
        status |= gr_poly_sub(t, t, u, R);
        status |= gr_poly_mul(res, S->psi + SLOT(S, j, m), t, R);
    }
    else
    {
        /* Psi_{m+2} Psi_m^3 and Psi_{m-1} Psi_{m+1}^3, one of them times W^2 */
        status |= gr_poly_mul(t, S->psi + SLOT(S, j, m + 2),
                    _divpoly_cb(S, j, m, R, &status), R);
        status |= gr_poly_mul(u, S->psi + SLOT(S, j, m - 1),
                    _divpoly_cb(S, j, m + 1, R, &status), R);

        if (m % 2)
        {
            status |= gr_poly_mul(res, W2, u, R);
            status |= gr_poly_sub(res, t, res, R);
        }
        else
        {
            status |= gr_poly_mul(res, W2, t, R);
            status |= gr_poly_sub(res, res, u, R);
        }
    }

    gr_poly_clear(t, R);
    gr_poly_clear(u, R);

    return status;
}

/* drop everything held at level j */
static void
_divpoly_free_level(divpoly_state_struct * S, slong j, gr_ctx_t R)
{
    slong s;

    for (s = S->off[j]; s < S->off[j + 1]; s++)
    {
        gr_poly_clear(S->psi + s, R);
        gr_poly_clear(S->sq + s, R);
        gr_poly_clear(S->cb + s, R);
        gr_poly_init(S->psi + s, R);
        gr_poly_init(S->sq + s, R);
        gr_poly_init(S->cb + s, R);
        S->state[s] = 0;
    }
}

int
gr_ec_ctx_division_poly(gr_poly_t res, ulong n, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    divpoly_state_struct S[1];
    gr_poly_t W, W2;
    slong j, s;
    ulong k;
    int status = GR_SUCCESS;

    if (n <= 4)
        return _divpoly_base(res, n, ctx);

    /* the result has degree about n^2/2 */
    if (n > (ulong) (WORD_MAX / n))
        return GR_UNABLE;

    /* index intervals, top-down, until everything is a base case */
    S->nlev = 0;
    S->lo[0] = S->hi[0] = n;

    while (S->hi[S->nlev] > 4)
    {
        ulong lo = S->lo[S->nlev], hi = S->hi[S->nlev];

        S->nlev++;
        S->lo[S->nlev] = (lo / 2 > 2) ? lo / 2 - 2 : 0;
        S->hi[S->nlev] = hi / 2 + 2;
    }

    S->nlev++;

    S->off[0] = 0;
    for (j = 0; j < S->nlev; j++)
        S->off[j + 1] = S->off[j] + (slong) (S->hi[j] - S->lo[j] + 1);

    S->total = S->off[S->nlev];

    S->psi = flint_malloc(S->total * sizeof(gr_poly_struct));
    S->sq = flint_malloc(S->total * sizeof(gr_poly_struct));
    S->cb = flint_malloc(S->total * sizeof(gr_poly_struct));
    S->state = flint_calloc(S->total, sizeof(unsigned char));

    for (s = 0; s < S->total; s++)
    {
        gr_poly_init(S->psi + s, R);
        gr_poly_init(S->sq + s, R);
        gr_poly_init(S->cb + s, R);
    }

    /* mark exactly the indices each level has to supply */
    S->state[SLOT(S, 0, n)] |= NEEDED;

    for (j = 0; j + 1 < S->nlev; j++)
    {
        for (k = S->lo[j]; k <= S->hi[j]; k++)
        {
            ulong m = k / 2, d, dlo;

            if (!(S->state[SLOT(S, j, k)] & NEEDED) || k <= 4)
                continue;

            dlo = (k % 2 == 0) ? m - 2 : m - 1;

            for (d = dlo; d <= m + 2; d++)
                S->state[SLOT(S, j + 1, d)] |= NEEDED;
        }
    }

    gr_poly_init(W, R);
    gr_poly_init(W2, R);
    status |= gr_ec_ctx_psi2_sqr(W, ctx);
    status |= gr_poly_mul(W2, W, W, R);

    /* evaluate bottom-up */
    for (j = S->nlev - 1; j >= 0 && status == GR_SUCCESS; j--)
    {
        for (k = S->lo[j]; k <= S->hi[j] && status == GR_SUCCESS; k++)
        {
            slong slot = SLOT(S, j, k);

            if (!(S->state[slot] & NEEDED))
                continue;

            if (k <= 4)
                status |= _divpoly_base(S->psi + slot, k, ctx);
            else
                status |= _divpoly_step(S->psi + slot, S, j + 1, k, W2, ctx);

            S->state[slot] |= HAVE_PSI;
        }

        /* the level below is no longer reachable */
        if (j + 1 < S->nlev)
            _divpoly_free_level(S, j + 1, R);
    }

    if (status == GR_SUCCESS)
        gr_poly_swap(res, S->psi + SLOT(S, 0, n), R);

    for (s = 0; s < S->total; s++)
    {
        gr_poly_clear(S->psi + s, R);
        gr_poly_clear(S->sq + s, R);
        gr_poly_clear(S->cb + s, R);
    }

    flint_free(S->psi);
    flint_free(S->sq);
    flint_free(S->cb);
    flint_free(S->state);

    gr_poly_clear(W, R);
    gr_poly_clear(W2, R);

    return status;
}

int
gr_ec_ctx_division_poly_vec(gr_poly_struct * res, slong len, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_t W, W2, t, u, v;
    slong k;
    int status = GR_SUCCESS;

    if (len <= 0)
        return GR_SUCCESS;

    for (k = 0; k < len && k <= 4; k++)
        status |= _divpoly_base(res + k, (ulong) k, ctx);

    if (len <= 5 || status != GR_SUCCESS)
        return status;

    gr_poly_init(W, R);
    gr_poly_init(W2, R);
    gr_poly_init(t, R);
    gr_poly_init(u, R);
    gr_poly_init(v, R);

    status |= gr_ec_ctx_psi2_sqr(W, ctx);
    status |= gr_poly_mul(W2, W, W, R);

    for (k = 5; k < len && status == GR_SUCCESS; k++)
    {
        slong m = k / 2;

        if (k % 2 == 0)
        {
            status |= gr_poly_mul(v, res + (m - 1), res + (m - 1), R);
            status |= gr_poly_mul(t, res + (m + 2), v, R);
            status |= gr_poly_mul(v, res + (m + 1), res + (m + 1), R);
            status |= gr_poly_mul(u, res + (m - 2), v, R);
            status |= gr_poly_sub(t, t, u, R);
            status |= gr_poly_mul(res + k, res + m, t, R);
        }
        else
        {
            status |= gr_poly_mul(v, res + m, res + m, R);
            status |= gr_poly_mul(v, v, res + m, R);
            status |= gr_poly_mul(t, res + (m + 2), v, R);

            status |= gr_poly_mul(v, res + (m + 1), res + (m + 1), R);
            status |= gr_poly_mul(v, v, res + (m + 1), R);
            status |= gr_poly_mul(u, res + (m - 1), v, R);

            if (m % 2)
            {
                status |= gr_poly_mul(v, W2, u, R);
                status |= gr_poly_sub(res + k, t, v, R);
            }
            else
            {
                status |= gr_poly_mul(v, W2, t, R);
                status |= gr_poly_sub(res + k, v, u, R);
            }
        }
    }

    gr_poly_clear(W, R);
    gr_poly_clear(W2, R);
    gr_poly_clear(t, R);
    gr_poly_clear(u, R);
    gr_poly_clear(v, R);

    return status;
}
