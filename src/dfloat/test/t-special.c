/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "gr.h"
#include "dfloat.h"

/* specific semantics: exactness, division by zero, infinite radii,
   overflow and underflow */

#define CHECK(cond) \
    do { if (!(cond)) { flint_printf("FAIL: %s (line %d)\n", #cond, __LINE__); flint_abort(); } } while (0)

TEST_FUNCTION_START(special, state)
{
    d2b_t x, y, z;
    d4b_t a, b, c;
    d1b_t p, q, r;
    gr_ctx_t ctx;
    int status;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    /* exactness */
    d2b_one(x);
    d2b_add(z, x, x);
    CHECK(z->d[0] == 2.0 && z->d[1] == 0.0 && z->rad == 0.0);
    d2b_set_d(y, 3.0);
    d2b_mul(z, y, y);
    CHECK(z->d[0] == 9.0 && z->rad == 0.0);
    d2b_set_dd(x, 0x1p60, 1.0);          /* 2^60 + 1 */
    d2b_mul(z, x, y);
    CHECK(z->d[0] == 3.0 * 0x1p60 && z->d[1] == 3.0 && z->rad == 0.0);
    d2b_sqr(z, x);                       /* needs 3 terms: inexact */
    CHECK(z->rad > 0.0 && z->rad <= 2.0);
    d2b_set_d(y, 1.0);
    d2b_set_d(x, 3.0);
    status = d2b_div(z, y, x);           /* 1/3 is inexact */
    CHECK(status == GR_SUCCESS && z->rad > 0.0 && z->rad < 0x1p-100);
    d2b_set_d(y, 6.0);
    status = d2b_div(z, y, x);           /* 6/3 = 2 exactly */
    CHECK(status == GR_SUCCESS && z->d[0] == 2.0 && z->rad == 0.0);
    d2b_set_d(y, 4.0);
    status = d2b_sqrt(z, y);
    CHECK(status == GR_SUCCESS && z->d[0] == 2.0 && z->rad == 0.0);
    d2b_set_d(y, 2.0);
    status = d2b_sqrt(z, y);
    CHECK(status == GR_SUCCESS && z->rad > 0.0 && z->rad < 0x1p-100);

    /* division by zero and by balls containing zero */
    d2b_zero(x);
    d2b_one(y);
    CHECK(d2b_div(z, y, x) == GR_DOMAIN);
    CHECK(d2b_inv(z, x) == GR_DOMAIN);
    d2b_set_d_rad(x, 0.0, 0.1);
    CHECK(d2b_div(z, y, x) == GR_UNABLE);
    d2b_set_d_rad(x, 0.05, 0.1);
    CHECK(d2b_div(z, y, x) == GR_UNABLE);
    d2b_set_d_rad(x, -1.0, 0.1);
    CHECK(d2b_div(z, y, x) == GR_SUCCESS);
    CHECK(d2b_contains_d(z, -1.0 / 1.1) && d2b_contains_d(z, -1.0 / 0.9));
    /* zero divided by anything nonzero is exactly zero */
    d2b_zero(y);
    CHECK(d2b_div(z, y, x) == GR_SUCCESS && d2b_is_zero(z));
    /* sqrt of negative and mixed-sign balls */
    d2b_set_d(x, -1.0);
    CHECK(d2b_sqrt(z, x) == GR_DOMAIN);
    d2b_set_d_rad(x, 0.0, 1.0);
    CHECK(d2b_sqrt(z, x) == GR_UNABLE);
    d2b_zero(x);
    CHECK(d2b_sqrt(z, x) == GR_SUCCESS && d2b_is_zero(z));

    /* infinite radius */
    d2b_indeterminate(x);
    d2b_one(y);
    d2b_add(z, x, y);
    CHECK(z->rad == D_INF && d2b_is_finite(z) == 0);
    d2b_mul(z, x, y);
    CHECK(z->rad == D_INF);
    d2b_zero(y);
    d2b_mul(z, x, y);                    /* exact zero times anything */
    CHECK(d2b_is_zero(z));
    d2b_mul(z, y, x);
    CHECK(d2b_is_zero(z));
    CHECK(d2b_div(z, y, x) == GR_UNABLE);
    CHECK(d2b_div(z, x, x) == GR_UNABLE);
    d2b_set_d_rad(y, 1.0, 0.5);
    status = d2b_div(z, x, y);
    CHECK(status == GR_SUCCESS && z->rad == D_INF);
    CHECK(d2b_sqrt(z, x) == GR_UNABLE);

    /* overflow: the whole real line; underflow: a tiny radius */
    d2b_set_d(x, DBL_MAX);
    d2b_add(z, x, x);
    CHECK(z->rad == D_INF && z->d[0] == 0.0);
    d2b_mul(z, x, x);
    CHECK(z->rad == D_INF);
    d2b_set_d(x, 0x1p-600);
    d2b_mul(z, x, x);
    CHECK(z->d[0] == 0.0 && z->rad > 0.0 && z->rad < 0x1p-1000);
    CHECK(d2b_contains_d(z, 0.0));
    d2b_set_d(x, 0x1p-540);
    d2b_mul(z, x, x);                    /* 2^-1080 is not representable */
    CHECK(z->rad > 0.0 && z->rad < 0x1p-1060);
    d2b_set_d(x, 0x1p-500);
    d2b_mul(z, x, x);                    /* 2^-1000 is */
    CHECK(z->d[0] == 0x1p-1000 && z->rad == 0.0);
    d2b_set_dd(x, 1.0, 0x1p-1070);       /* a legitimate but extreme tail */
    d2b_set_d(y, 3.0);
    d2b_mul(z, x, y);
    CHECK(d2b_contains_d(z, 3.0) && z->rad > 0.0 && z->rad < 0x1p-1000);
    d2b_set_d(x, 0x1p1000);
    d2b_set_d(y, 0x1p-1000);
    status = d2b_div(z, y, x);
    CHECK(status == GR_SUCCESS && z->d[0] == 0.0 && z->rad > 0.0 && z->rad < 0x1p-1060);
    status = d2b_div(z, x, y);
    CHECK(status == GR_SUCCESS && z->rad == D_INF);

    /* d1b and d4b sanity */
    d1b_set_d(p, 0.1);
    d1b_set_d(q, 0.2);
    d1b_add(r, p, q);
    {
        /* the sum of the doubles nearest to 0.1 and 0.2 (written with
           variables: with x87 excess precision the constant expression
           0.1 + 0.2 may be evaluated as 0.1L + 0.2L) */
        double t1 = 0.1, t2 = 0.2;
        CHECK(r->rad > 0.0 && d1b_contains_d(r, t1 + t2) && r->rad < 0x1p-52);
    }
    d1b_mul(r, p, q);
    CHECK(r->rad > 0.0 && r->rad < 0x1p-56);
    d4b_set_d(a, 0.1);
    d4b_set_d(b, 0.2);
    d4b_add(c, a, b);
    CHECK(c->rad == 0.0);                /* representable in two doubles */
    d4b_mul(c, a, b);
    CHECK(c->rad == 0.0);                /* 106 bits fit in four doubles */
    d4b_sqr(c, c);
    CHECK(c->rad == 0.0);                /* 212 bits, still fits */
    d4b_mul(c, c, a);
    CHECK(c->rad > 0.0 && c->rad < 0x1p-225);

    /* gr layer: division by zero */
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, 3, 1));
    {
        gr_ptr u, v;
        GR_TMP_INIT2(u, v, ctx);
        GR_MUST_SUCCEED(gr_one(u, ctx));
        GR_MUST_SUCCEED(gr_zero(v, ctx));
        CHECK(gr_div(v, u, v, ctx) == GR_DOMAIN);
        CHECK(gr_inv(v, v, ctx) == GR_DOMAIN);
        CHECK(gr_is_invertible(v, ctx) == T_FALSE);
        CHECK(gr_is_zero(v, ctx) == T_TRUE);
        CHECK(gr_equal(u, v, ctx) == T_FALSE);
        d3b_set_d_rad(v, 0.0, 0.1);
        CHECK(gr_div(v, u, v, ctx) == GR_UNABLE);
        CHECK(gr_is_zero(v, ctx) == T_UNKNOWN);
        CHECK(gr_is_invertible(v, ctx) == T_UNKNOWN);
        GR_TMP_CLEAR2(u, v, ctx);
    }
    gr_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
