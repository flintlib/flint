/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "arb.h"
#include "gr.h"

/* random ball with moderate exponents; stores a random rational point of
   the ball in q */
static void
_random_ball_point(decball_t x, fmpq_t q, flint_rand_t state, gr_ctx_t ctx)
{
    fmpq_t m, r, t;

    fmpq_init(m);
    fmpq_init(r);
    fmpq_init(t);

    do
    {
        GR_MUST_SUCCEED(decfloat_randtest_special(&x->mid, state, ctx));
    }
    while (!DECFLOAT_IS_FINITE(&x->mid));

    if (!DECFLOAT_IS_SPECIAL(&x->mid))
        fmpz_set_si(&x->mid.exp, (slong) n_randint(state, 21) - 10);

    switch (n_randint(state, 4))
    {
        case 0:
            _decmag_zero(&x->rad, ctx);
            break;
        case 1:
            {
                slong E;
                if (!DECFLOAT_IS_SPECIAL(&x->mid) && decfloat_get_sci_exp_si(&E, &x->mid, ctx))
                    _decmag_set_ui_10exp_si(&x->rad, 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1),
                        E - DECIMAL_CTX_PREC(ctx) - DECIMAL_CTX_RAD_PREC(ctx) + n_randint(state, 4), ctx);
                else
                    _decmag_zero(&x->rad, ctx);
            }
            break;
        default:
            _decmag_set_ui_10exp_si(&x->rad, 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1),
                (slong) n_randint(state, 41) - 30, ctx);
            break;
    }

    GR_MUST_SUCCEED(decfloat_get_fmpq(m, &x->mid, ctx));
    _decmag_get_fmpq(r, &x->rad, ctx);

    /* q = m + r * s with s in [-1, 1] */
    switch (n_randint(state, 6))
    {
        case 0: fmpq_set_si(t, 1, 1); break;
        case 1: fmpq_set_si(t, -1, 1); break;
        case 2: fmpq_zero(t); break;
        default:
            {
                fmpz_t a, b;
                fmpz_init(a);
                fmpz_init(b);
                fmpz_randtest(a, state, 30);
                fmpz_randtest_not_zero(b, state, 30);
                fmpz_abs(b, b);
                fmpz_fdiv_r(a, a, b);   /* 0 <= a < b */
                fmpq_set_fmpz_frac(t, a, b);
                if (n_randint(state, 2)) fmpq_neg(t, t);
                fmpz_clear(a);
                fmpz_clear(b);
            }
            break;
    }

    fmpq_mul(t, t, r);
    fmpq_add(q, m, t);

    fmpq_clear(m);
    fmpq_clear(r);
    fmpq_clear(t);
}

TEST_FUNCTION_START(decball, state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ctx_init_decball_randtest(ctx, state, 40);
        gr_test_ring(ctx, 30, 0);
        gr_ctx_clear(ctx);
    }

    {
        gr_ctx_t ctx;
        gr_ctx_init_decball(ctx, 10, 0);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
        gr_ctx_init_decball(ctx, 30, DECIMAL_SLOPPY_RADIUS);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
    }

    /* containment */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decball_t X, Y, Z;
        fmpq_t x, y, z;
        int op, status, aliasing;
        slong prec;

        gr_ctx_init_decball_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        prec = DECIMAL_CTX_PREC(ctx);

        decball_init(X, ctx);
        decball_init(Y, ctx);
        decball_init(Z, ctx);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        _random_ball_point(X, x, state, ctx);
        _random_ball_point(Y, y, state, ctx);

        if (!_decball_contains_fmpq(X, x, ctx) || !_decball_contains_fmpq(Y, y, ctx))
        {
            flint_printf("FAIL: contains_fmpq\n");
            flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
            flint_printf("Y = %{gr}, y = %{fmpq}\n", Y, ctx, y);
            flint_abort();
        }

        op = n_randint(state, 12);
        aliasing = n_randint(state, 3);

        if (aliasing == 1)
        {
            GR_MUST_SUCCEED(decball_set(Z, X, ctx));
            /* rounding may have enlarged the ball; keep it consistent */
            GR_MUST_SUCCEED(decball_set(X, Z, ctx));
        }
        else if (aliasing == 2)
        {
            GR_MUST_SUCCEED(decball_set(Z, Y, ctx));
            GR_MUST_SUCCEED(decball_set(Y, Z, ctx));
        }

        {
            const decball_struct * XA = (aliasing == 1) ? Z : X;
            const decball_struct * YA = (aliasing == 2) ? Z : Y;

            switch (op)
            {
                case 0: status = decball_add(Z, XA, YA, ctx); fmpq_add(z, x, y); break;
                case 1: status = decball_sub(Z, XA, YA, ctx); fmpq_sub(z, x, y); break;
                case 2: status = decball_mul(Z, XA, YA, ctx); fmpq_mul(z, x, y); break;
                case 3:
                    status = decball_div(Z, XA, YA, ctx);
                    if (fmpq_is_zero(y))
                    {
                        if (status == GR_SUCCESS)
                        {
                            flint_printf("FAIL: division by a ball containing zero succeeded\n");
                            flint_printf("Y = %{gr}\n", Y, ctx);
                            flint_abort();
                        }
                        fmpq_zero(z);
                    }
                    else
                        fmpq_div(z, x, y);
                    break;
                case 4: status = decball_neg(Z, XA, ctx); fmpq_neg(z, x); break;
                case 5: status = decball_abs(Z, XA, ctx); fmpq_abs(z, x); break;
                case 6: status = decball_sqr(Z, XA, ctx); fmpq_mul(z, x, x); break;
                case 7:
                    status = decball_inv(Z, XA, ctx);
                    if (fmpq_is_zero(x))
                    {
                        if (status == GR_SUCCESS)
                        {
                            flint_printf("FAIL: inverse of a ball containing zero succeeded\n");
                            flint_abort();
                        }
                        fmpq_zero(z);
                    }
                    else
                        fmpq_inv(z, x);
                    break;
                case 8:
                    /* sqrt: test via z^2 containment below */
                    status = decball_sqrt(Z, XA, ctx);
                    fmpq_set(z, x);
                    if (fmpq_sgn(x) < 0 && status == GR_SUCCESS)
                    {
                        flint_printf("FAIL: sqrt of a ball containing negative numbers succeeded\n");
                        flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                        flint_abort();
                    }
                    break;
                case 9: status = decball_set(Z, XA, ctx); fmpq_set(z, x); break;
                case 10:
                    {
                        slong e = (slong) n_randint(state, 41) - 20;
                        fmpq_t p;
                        fmpq_init(p);
                        fmpq_one(p);
                        if (e >= 0) fmpz_ui_pow_ui(fmpq_numref(p), 10, e);
                        else fmpz_ui_pow_ui(fmpq_denref(p), 10, -e);
                        status = decball_mul_10exp_si(Z, XA, e, ctx);
                        fmpq_mul(z, x, p);
                        fmpq_clear(p);
                    }
                    break;
                default:
                    {
                        /* floor/ceil/trunc/nint */
                        fmpz_t n;
                        fmpz_init(n);
                        switch (n_randint(state, 4))
                        {
                            case 0: status = decball_floor(Z, XA, ctx); fmpz_fdiv_q(n, fmpq_numref(x), fmpq_denref(x)); break;
                            case 1: status = decball_ceil(Z, XA, ctx); fmpz_cdiv_q(n, fmpq_numref(x), fmpq_denref(x)); break;
                            case 2: status = decball_trunc(Z, XA, ctx); fmpz_tdiv_q(n, fmpq_numref(x), fmpq_denref(x)); break;
                            default:
                                status = decball_nint(Z, XA, ctx);
                                {
                                    /* nearest with ties to even */
                                    fmpq_t h;
                                    fmpz_t r2;
                                    fmpq_init(h);
                                    fmpz_init(r2);
                                    fmpz_fdiv_qr(n, r2, fmpq_numref(x), fmpq_denref(x));
                                    fmpq_set_fmpz_frac(h, r2, fmpq_denref(x));   /* fractional part in [0, 1) */
                                    fmpq_mul_2exp(h, h, 1);
                                    if (fmpq_cmp_si(h, 1) > 0 || (fmpq_cmp_si(h, 1) == 0 && fmpz_is_odd(n)))
                                        fmpz_add_ui(n, n, 1);
                                    fmpq_clear(h);
                                    fmpz_clear(r2);
                                }
                                break;
                        }
                        fmpq_set_fmpz(z, n);
                        fmpz_clear(n);
                    }
                    break;
            }
        }

        if (status == GR_SUCCESS)
        {
            int ok;

            if (op == 8)
            {
                /* z = x >= 0: Z must contain sqrt(x): check via the interval [lo, hi]
                   of Z: lo^2 <= x <= hi^2 (with lo >= 0) */
                fmpq_t m, r, lo, hi;
                fmpq_init(m); fmpq_init(r); fmpq_init(lo); fmpq_init(hi);
                GR_MUST_SUCCEED(decfloat_get_fmpq(m, &Z->mid, ctx));
                _decmag_get_fmpq(r, &Z->rad, ctx);
                fmpq_sub(lo, m, r);
                fmpq_add(hi, m, r);
                if (fmpq_sgn(lo) < 0) fmpq_zero(lo);
                fmpq_mul(lo, lo, lo);
                fmpq_mul(hi, hi, hi);
                ok = (fmpq_cmp(lo, z) <= 0 && fmpq_cmp(z, hi) <= 0);
                fmpq_clear(m); fmpq_clear(r); fmpq_clear(lo); fmpq_clear(hi);
            }
            else
            {
                ok = _decball_contains_fmpq(Z, z, ctx);
            }

            if (!ok)
            {
                flint_printf("FAIL: containment, op = %d, aliasing = %d\n", op, aliasing);
                gr_ctx_println(ctx);
                flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                flint_printf("Y = %{gr}, y = %{fmpq}\n", Y, ctx, y);
                flint_printf("Z = %{gr}, z = %{fmpq}\n", Z, ctx, z);
                flint_abort();
            }

            /* radius sanity: for exact inputs, the result is at most about 1 ulp wide */
            if (op <= 2 && DECMAG_IS_ZERO(&X->rad) && DECMAG_IS_ZERO(&Y->rad) && !DECFLOAT_IS_SPECIAL(&Z->mid))
            {
                decmag_t u;
                _decmag_init(u, ctx);
                _decmag_set_ulp(u, &Z->mid, prec, ctx);
                if (_decmag_cmp(&Z->rad, u, ctx) > 0)
                {
                    flint_printf("FAIL: radius too large, op = %d\n", op);
                    gr_ctx_println(ctx);
                    flint_printf("X = %{gr}\n", X, ctx);
                    flint_printf("Y = %{gr}\n", Y, ctx);
                    flint_printf("Z = %{gr}\n", Z, ctx);
                    flint_abort();
                }
                _decmag_clear(u, ctx);
            }
        }
        else if (op <= 2 || op == 4 || op == 5 || op == 6 || op == 9 || op == 10)
        {
            flint_printf("FAIL: unexpected failure, op = %d\n", op);
            gr_ctx_println(ctx);
            flint_printf("X = %{gr}\n", X, ctx);
            flint_printf("Y = %{gr}\n", Y, ctx);
            flint_abort();
        }

        /* predicates */
        {
            truth_t t;
            int c;

            t = decball_equal(X, Y, ctx);
            if ((t == T_TRUE && !fmpq_equal(x, y)) || (t == T_FALSE && fmpq_equal(x, y)))
            {
                flint_printf("FAIL: equal\n");
                flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                flint_printf("Y = %{gr}, y = %{fmpq}\n", Y, ctx, y);
                flint_abort();
            }

            if (decball_cmp(&c, X, Y, ctx) == GR_SUCCESS && c != fmpq_cmp(x, y))
            {
                flint_printf("FAIL: cmp\n");
                flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                flint_printf("Y = %{gr}, y = %{fmpq}\n", Y, ctx, y);
                flint_abort();
            }

            t = decball_is_zero(X, ctx);
            if ((t == T_TRUE && !fmpq_is_zero(x)) || (t == T_FALSE && fmpq_is_zero(x)))
            {
                flint_printf("FAIL: is_zero\n");
                flint_abort();
            }

            if (_decball_is_positive(X, ctx) && fmpq_sgn(x) <= 0) { flint_printf("FAIL: is_positive\n"); flint_abort(); }
            if (_decball_is_negative(X, ctx) && fmpq_sgn(x) >= 0) { flint_printf("FAIL: is_negative\n"); flint_abort(); }
            if (_decball_is_nonnegative(X, ctx) && fmpq_sgn(x) < 0) { flint_printf("FAIL: is_nonnegative\n"); flint_abort(); }
            if (_decball_is_nonpositive(X, ctx) && fmpq_sgn(x) > 0) { flint_printf("FAIL: is_nonpositive\n"); flint_abort(); }

            if (!_decball_overlaps(X, X, ctx)) { flint_printf("FAIL: overlaps self\n"); flint_abort(); }
            if (!_decball_contains(X, X, ctx)) { flint_printf("FAIL: contains self\n"); flint_abort(); }

            if (fmpq_equal(x, y) && !_decball_overlaps(X, Y, ctx))
            {
                flint_printf("FAIL: overlaps\n");
                flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                flint_printf("Y = %{gr}, y = %{fmpq}\n", Y, ctx, y);
                flint_abort();
            }

            t = decball_is_integer(X, ctx);
            if ((t == T_TRUE && !fmpz_is_one(fmpq_denref(x))) || (t == T_FALSE && fmpz_is_one(fmpq_denref(x))))
            {
                flint_printf("FAIL: is_integer\n");
                flint_printf("X = %{gr}, x = %{fmpq}\n", X, ctx, x);
                flint_abort();
            }
        }

        decball_clear(X, ctx);
        decball_clear(Y, ctx);
        decball_clear(Z, ctx);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        gr_ctx_clear(ctx);
    }

    /* contains and overlaps against exact rational arithmetic, including
       tight cases where the radii are exactly the distance of the midpoints */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decball_t X, Y;
        fmpq_t xm, ym, xr, yr, d, t;
        int c1, c2, o1, o2;

        gr_ctx_init_decball_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decball_init(X, ctx);
        decball_init(Y, ctx);
        fmpq_init(xm); fmpq_init(ym); fmpq_init(xr); fmpq_init(yr); fmpq_init(d); fmpq_init(t);

        GR_MUST_SUCCEED(decball_randtest(X, state, ctx));
        GR_MUST_SUCCEED(decball_randtest(Y, state, ctx));

        /* keep the exponents moderate so that the rational oracle is cheap */
        if (fmpz_bits(&X->mid.exp) > 6 || fmpz_bits(&Y->mid.exp) > 6 ||
            fmpz_bits(&X->rad.exp) > 8 || fmpz_bits(&Y->rad.exp) > 8)
        {
            decball_clear(X, ctx); decball_clear(Y, ctx);
            fmpq_clear(xm); fmpq_clear(ym); fmpq_clear(xr); fmpq_clear(yr); fmpq_clear(d); fmpq_clear(t);
            gr_ctx_clear(ctx);
            continue;
        }

        switch (n_randint(state, 6))
        {
            case 0:
                /* Y = X shifted, same radius */
                decball_set(Y, X, ctx);
                GR_MUST_SUCCEED(decball_add_si(Y, Y, (slong) n_randint(state, 3) - 1, ctx));
                break;
            case 1:
                /* Y = point at the edge of X (up to rounding) */
                {
                    decfloat_t r;
                    decfloat_init(r, ctx);
                    GR_MUST_SUCCEED(_decmag_get_decfloat(r, &X->rad, ctx));
                    if (n_randint(state, 2)) r->m.size = -r->m.size;
                    decimal_ctx_set_prec(ctx, DECIMAL_PREC_EXACT);
                    GR_MUST_SUCCEED(decfloat_add(&Y->mid, &X->mid, r, ctx));
                    _decmag_zero(&Y->rad, ctx);
                    if (n_randint(state, 2))
                        _decmag_set(&Y->rad, &X->rad, ctx);
                    decfloat_clear(r, ctx);
                }
                break;
            case 2:
                /* radius of X = distance to Y plus radius of Y (tight) */
                {
                    decfloat_t r;
                    decfloat_init(r, ctx);
                    decimal_ctx_set_prec(ctx, DECIMAL_PREC_EXACT);
                    GR_MUST_SUCCEED(decfloat_sub(r, &X->mid, &Y->mid, ctx));
                    _decmag_set_decfloat(&X->rad, r, ctx);
                    if (n_randint(state, 2))
                        _decmag_add(&X->rad, &X->rad, &Y->rad, ctx);
                    else
                        _decmag_add_lower(&X->rad, &X->rad, &Y->rad, ctx);
                    decfloat_clear(r, ctx);
                }
                break;
            default:
                break;
        }

        GR_MUST_SUCCEED(decfloat_get_fmpq(xm, &X->mid, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(ym, &Y->mid, ctx));
        _decmag_get_fmpq(xr, &X->rad, ctx);
        _decmag_get_fmpq(yr, &Y->rad, ctx);

        fmpq_sub(d, xm, ym);
        fmpq_abs(d, d);

        /* contains: d + yr <= xr */
        fmpq_add(t, d, yr);
        c1 = (fmpq_cmp(t, xr) <= 0);
        c2 = _decball_contains(X, Y, ctx);

        /* overlaps: d <= xr + yr */
        fmpq_add(t, xr, yr);
        o1 = (fmpq_cmp(d, t) <= 0);
        o2 = _decball_overlaps(X, Y, ctx);

        if (c1 != c2 || o1 != o2)
        {
            flint_printf("FAIL: contains / overlaps\n");
            gr_ctx_println(ctx);
            flint_printf("X = %{gr}\n", X, ctx);
            flint_printf("Y = %{gr}\n", Y, ctx);
            flint_printf("contains: exact %d, got %d\n", c1, c2);
            flint_printf("overlaps: exact %d, got %d\n", o1, o2);
            flint_abort();
        }

        decball_clear(X, ctx);
        decball_clear(Y, ctx);
        fmpq_clear(xm); fmpq_clear(ym); fmpq_clear(xr); fmpq_clear(yr); fmpq_clear(d); fmpq_clear(t);
        gr_ctx_clear(ctx);
    }

    /* endpoints, components, hulls and rounding of the radius; the operands
       may come from a context with a different radius precision */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, ctx0;
        decball_t X, Y, Z, W;
        fmpq_t xm, xr, a, b, c;
        slong rp2;

        gr_ctx_init_decball_randtest(ctx, state, 40);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        /* same limb radix, other precision, radius precision and flags */
        _gr_ctx_init_decimal(ctx0, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), 1 + n_randint(state, 40), DECIMAL_RND_DOWN,
            n_randint(state, 2) ? DECIMAL_SLOPPY_RADIUS : 0);
        decimal_ctx_set_rad_prec(ctx0, 1 + n_randint(state, DECMAG_MAX_PREC));

        decball_init(X, ctx);
        decball_init(Y, ctx);
        decball_init(Z, ctx);
        decball_init(W, ctx);
        fmpq_init(xm); fmpq_init(xr); fmpq_init(a); fmpq_init(b); fmpq_init(c);

        /* X and Y are created in ctx0 */
        _random_ball_point(X, a, state, ctx0);
        _random_ball_point(Y, b, state, ctx0);
        if (DECFLOAT_IS_SPECIAL(&X->mid))
            GR_MUST_SUCCEED(decball_one(X, ctx0));
        if (DECFLOAT_IS_SPECIAL(&Y->mid))
        {
            GR_MUST_SUCCEED(decball_one(Y, ctx0));
            fmpq_one(b);
        }

        GR_MUST_SUCCEED(decfloat_get_fmpq(xm, &X->mid, ctx));
        _decmag_get_fmpq(xr, &X->rad, ctx);

        /* lower <= xm - xr, upper >= xm + xr, both exact balls */
        GR_MUST_SUCCEED(decball_lower(Z, X, ctx));
        GR_MUST_SUCCEED(decball_upper(W, X, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(a, &Z->mid, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(b, &W->mid, ctx));
        fmpq_sub(c, xm, xr);
        if (!_decball_is_exact(Z, ctx) || !_decball_is_exact(W, ctx) || fmpq_cmp(a, c) > 0)
        {
            flint_printf("FAIL: lower\nX = %{gr}\nZ = %{gr}\n", X, ctx, Z, ctx);
            flint_abort();
        }
        fmpq_add(c, xm, xr);
        if (fmpq_cmp(b, c) < 0)
        {
            flint_printf("FAIL: upper\nX = %{gr}\nW = %{gr}\n", X, ctx, W, ctx);
            flint_abort();
        }
        /* the endpoints are within one ulp at the context precision:
           set_interval(lower, upper) must contain X and be at most
           slightly wider */
        GR_MUST_SUCCEED(decball_set_interval(Z, Z, W, ctx));
        if (!_decball_contains(Z, X, ctx))
        {
            flint_printf("FAIL: set_interval(lower, upper) contains\nX = %{gr}\nZ = %{gr}\n", X, ctx, Z, ctx);
            flint_abort();
        }

        /* abs bounds */
        GR_MUST_SUCCEED(decball_abs_lower(Z, X, ctx));
        GR_MUST_SUCCEED(decball_abs_upper(W, X, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(a, &Z->mid, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(b, &W->mid, ctx));
        fmpq_abs(c, xm);
        fmpq_add(c, c, xr);      /* max |x| */
        if (fmpq_cmp(b, c) < 0 || fmpq_sgn(a) < 0)
        {
            flint_printf("FAIL: abs_upper\nX = %{gr}\nW = %{gr}\n", X, ctx, W, ctx);
            flint_abort();
        }
        fmpq_abs(c, xm);
        fmpq_sub(c, c, xr);      /* min |x| (if positive) */
        if (fmpq_sgn(c) < 0)
            fmpq_zero(c);
        if (fmpq_cmp(a, c) > 0)
        {
            flint_printf("FAIL: abs_lower\nX = %{gr}\nZ = %{gr}\n", X, ctx, Z, ctx);
            flint_abort();
        }

        /* mid + shell = X, rad = radius */
        GR_MUST_SUCCEED(decball_mid(Z, X, ctx));
        GR_MUST_SUCCEED(decball_shell(W, X, ctx));
        GR_MUST_SUCCEED(decball_rad(Y, X, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(a, &Y->mid, ctx));
        GR_MUST_SUCCEED(decfloat_get_fmpq(b, &Z->mid, ctx));
        if (!_decball_is_exact(Z, ctx) || !_decball_is_exact(Y, ctx) || !fmpq_equal(a, xr) || !fmpq_equal(b, xm)
            || !DECFLOAT_IS_ZERO(&W->mid) || _decmag_cmp(&W->rad, &X->rad, ctx) < 0 || !DECMAG_IS_NORMALIZED(&W->rad, ctx))
        {
            flint_printf("FAIL: mid/rad/shell\nX = %{gr}\nZ = %{gr}\nY = %{gr}\nW = %{gr}\n", X, ctx, Z, ctx, Y, ctx, W, ctx);
            flint_abort();
        }
        GR_MUST_SUCCEED(decball_add(Z, Z, W, ctx));
        if (!_decball_contains(Z, X, ctx))
        {
            flint_printf("FAIL: mid + shell\nX = %{gr}\nZ = %{gr}\n", X, ctx, Z, ctx);
            flint_abort();
        }

        /* add_rad: contains X shifted by any point of Y */
        _random_ball_point(Y, b, state, ctx0);
        if (DECFLOAT_IS_SPECIAL(&Y->mid))
        {
            GR_MUST_SUCCEED(decball_one(Y, ctx0));
            fmpq_one(b);
        }
        GR_MUST_SUCCEED(decball_add_rad(Z, X, Y, ctx));
        fmpq_add(c, xm, b);
        fmpq_sub(a, xm, b);
        if (!_decball_contains_fmpq(Z, c, ctx) || !_decball_contains_fmpq(Z, a, ctx) || !_decball_contains(Z, X, ctx))
        {
            flint_printf("FAIL: add_rad\nX = %{gr}\nY = %{gr}\nZ = %{gr}\n", X, ctx, Y, ctx, Z, ctx);
            flint_abort();
        }

        /* set_interval: hull of two balls (in either order) */
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(decball_set_interval(Z, X, Y, ctx));
        else
            GR_MUST_SUCCEED(decball_set_interval(Z, Y, X, ctx));
        if (!_decball_contains(Z, X, ctx) || !_decball_contains(Z, Y, ctx))
        {
            flint_printf("FAIL: set_interval\nX = %{gr}\nY = %{gr}\nZ = %{gr}\n", X, ctx, Y, ctx, Z, ctx);
            flint_abort();
        }

        /* set_round2: contains X, radius rounded to rp2 digits */
        rp2 = 1 + n_randint(state, DECMAG_MAX_PREC);
        GR_MUST_SUCCEED(decball_set_round2(Z, X, 1 + n_randint(state, 30), rp2, ctx));
        if (!_decball_contains(Z, X, ctx) || (!DECMAG_IS_SPECIAL(&Z->rad) && _decmag_digits(&Z->rad) != rp2))
        {
            flint_printf("FAIL: set_round2\nX = %{gr}\nZ = %{gr}\nrp2 = %wd\n", X, ctx, Z, ctx, rp2);
            flint_abort();
        }

        /* arithmetic with operands of another radius precision: the results
           are normalized to this context and contain the exact results */
        GR_MUST_SUCCEED(decball_set_interval_mid_rad(W, X, Y, ctx));
        GR_MUST_SUCCEED(decball_add(Z, X, Y, ctx));
        fmpq_add(c, xm, b);
        if (!_decball_contains_fmpq(Z, c, ctx) || !DECMAG_IS_NORMALIZED(&Z->rad, ctx) || !DECMAG_IS_NORMALIZED(&W->rad, ctx))
        {
            flint_printf("FAIL: add (mixed radius precision)\nX = %{gr}\nY = %{gr}\nZ = %{gr}\nW = %{gr}\n", X, ctx, Y, ctx, Z, ctx, W, ctx);
            flint_abort();
        }
        GR_MUST_SUCCEED(decball_mul(Z, X, Y, ctx));
        fmpq_mul(c, xm, b);
        if (!_decball_contains_fmpq(Z, c, ctx) || !DECMAG_IS_NORMALIZED(&Z->rad, ctx))
        {
            flint_printf("FAIL: mul (mixed radius precision)\nX = %{gr}\nY = %{gr}\nZ = %{gr}\n", X, ctx, Y, ctx, Z, ctx);
            flint_abort();
        }
        GR_MUST_SUCCEED(decball_set(Z, X, ctx));
        if (!_decball_contains_fmpq(Z, xm, ctx) || !DECMAG_IS_NORMALIZED(&Z->rad, ctx) || !_decball_contains(Z, X, ctx))
        {
            flint_printf("FAIL: set (mixed radius precision)\nX = %{gr}\nZ = %{gr}\n", X, ctx, Z, ctx);
            flint_abort();
        }

        decball_clear(X, ctx);
        decball_clear(Y, ctx);
        decball_clear(Z, ctx);
        decball_clear(W, ctx);
        fmpq_clear(xm); fmpq_clear(xr); fmpq_clear(a); fmpq_clear(b); fmpq_clear(c);
        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx0);
    }

    TEST_FUNCTION_END(state);
}
