/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "padic_radix.h"

/*
    The fused methods (addmul / submul and their scalar versions, the scalar
    multiplications, the vector-scalar operations) must give bit-identical
    results to the obvious compositions of the core operations, and pow must
    agree with one rounding of the exact power (with the error term
    (n - 1) v + N + v_p(n)).
*/

static int
_same(const padic_radix_t a, const padic_radix_t b, gr_ctx_t ctx)
{
    return a->v == b->v && a->N == b->N
        && radix_integer_equal(&a->u, &b->u, PADIC_RADIX_CTX_RADIX(ctx));
}

static void
_randelem(padic_radix_t x, flint_rand_t state, gr_ctx_t ctx)
{
    slong e = radix_limb_exponent(PADIC_RADIX_CTX_RADIX(ctx));

    padic_radix_randtest(x, state, ctx);

    /* valuations spread over a few limbs, and low relative precision */
    if (x->u.size != 0 && n_randint(state, 2))
    {
        slong rel = (x->N == PADIC_RADIX_EXACT) ? 0 : x->N - x->v;
        x->v = (slong) n_randint(state, 6 * e + 1) - 3 * e;
        if (x->N != PADIC_RADIX_EXACT)
            x->N = x->v + rel;
        if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
            padic_radix_zero(x, ctx);
    }

    if (x->u.size != 0 && n_randint(state, 3) == 0
        && (PADIC_RADIX_CTX_PREC_ABS(ctx) != PADIC_RADIX_PREC_INF
            || PADIC_RADIX_CTX_PREC_REL(ctx) != PADIC_RADIX_PREC_INF))
    {
        x->N = x->v + 1 + n_randint(state, 4 * e);
        if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
            padic_radix_zero(x, ctx);
    }
}

static slong
_randscalar(flint_rand_t state, ulong p)
{
    slong c;
    switch (n_randint(state, 4))
    {
        case 0:  c = (slong) n_randint(state, 21) - 10; break;
        case 1:  c = (slong) n_randtest(state); break;
        case 2:  c = n_randint(state, 2) ? WORD_MIN : WORD_MAX; break;
        default:
        {
            ulong t = 1 + n_randint(state, 20);
            while (n_randint(state, 3) && !n_mul_checked(&t, t, p))
                ;
            c = (slong) (t >> 1);
        }
    }
    return c;
}

#define CHECK(op, r1, r2) \
    do { \
        if (s1 != s2 || (s1 == GR_SUCCESS && !_same(r1, r2, ctx))) \
        { \
            flint_printf("FAIL: %s (s1 = %d, s2 = %d)\n", op, s1, s2); \
            gr_ctx_println(ctx); \
            flint_printf("x = "); gr_println(x, ctx); \
            flint_printf("y = "); gr_println(y, ctx); \
            flint_printf("r1 = "); gr_println(r1, ctx); \
            flint_printf("r2 = "); gr_println(r2, ctx); \
            flint_abort(); \
        } \
    } while (0)

TEST_FUNCTION_START(padic_radix_ops, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_t x, y, r0, r1, r2, t;
        fmpz_t cz;
        ulong p;
        slong c, i;
        int s1, s2, sub;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 200));
        p = GR_PADIC_RADIX_CTX(ctx)->p;

        padic_radix_init(x, ctx);
        padic_radix_init(y, ctx);
        padic_radix_init(r0, ctx);
        padic_radix_init(r1, ctx);
        padic_radix_init(r2, ctx);
        padic_radix_init(t, ctx);
        fmpz_init(cz);

        _randelem(x, state, ctx);
        _randelem(y, state, ctx);
        _randelem(r0, state, ctx);

        /* make cancellation likely */
        if (n_randint(state, 3) == 0 && padic_radix_mul(t, x, y, ctx) == GR_SUCCESS)
            GR_IGNORE(padic_radix_neg(r0, t, ctx));

        sub = n_randint(state, 2);

        /* addmul / submul vs mul + add */
        GR_IGNORE(padic_radix_set(r1, r0, ctx));
        s1 = padic_radix_mul(t, x, y, ctx);
        if (s1 == GR_SUCCESS)
            s1 = sub ? padic_radix_sub(r1, r1, t, ctx) : padic_radix_add(r1, r1, t, ctx);
        GR_IGNORE(padic_radix_set(r2, r0, ctx));
        s2 = sub ? padic_radix_submul(r2, x, y, ctx) : padic_radix_addmul(r2, x, y, ctx);
        CHECK("addmul / submul", r1, r2);

        /* aliasing res = x */
        GR_IGNORE(padic_radix_set(r1, x, ctx));
        s1 = padic_radix_mul(t, x, y, ctx);
        if (s1 == GR_SUCCESS)
            s1 = padic_radix_add(r1, r1, t, ctx);
        GR_IGNORE(padic_radix_set(r2, x, ctx));
        s2 = padic_radix_addmul(r2, r2, y, ctx);
        CHECK("addmul aliased", r1, r2);

        /* scalars */
        c = _randscalar(state, p);

        fmpz_set_si(cz, c);
        s1 = padic_radix_set_fmpz(r1, cz, ctx);
        s2 = padic_radix_set_si(r2, c, ctx);
        CHECK("set_si", r1, r2);

        fmpz_set_ui(cz, (ulong) c);
        s1 = padic_radix_set_fmpz(r1, cz, ctx);
        s2 = padic_radix_set_ui(r2, (ulong) c, ctx);
        CHECK("set_ui", r1, r2);

        fmpz_set_si(cz, c);
        s1 = padic_radix_set_fmpz(t, cz, ctx);
        if (s1 == GR_SUCCESS)
            s1 = padic_radix_mul(r1, x, t, ctx);
        s2 = padic_radix_mul_si(r2, x, c, ctx);
        CHECK("mul_si", r1, r2);

        s1 = padic_radix_set_fmpz(t, cz, ctx);
        GR_IGNORE(padic_radix_set(r1, r0, ctx));
        if (s1 == GR_SUCCESS)
            s1 = padic_radix_addmul(r1, x, t, ctx);
        GR_IGNORE(padic_radix_set(r2, r0, ctx));
        s2 = padic_radix_addmul_si(r2, x, c, ctx);
        CHECK("addmul_si", r1, r2);

        fmpz_set_ui(cz, (ulong) c);
        s1 = padic_radix_set_fmpz(t, cz, ctx);
        GR_IGNORE(padic_radix_set(r1, r0, ctx));
        if (s1 == GR_SUCCESS)
            s1 = padic_radix_submul(r1, x, t, ctx);
        GR_IGNORE(padic_radix_set(r2, r0, ctx));
        s2 = padic_radix_submul_ui(r2, x, (ulong) c, ctx);
        CHECK("submul_ui", r1, r2);

        fmpz_randtest(cz, state, 300);
        s1 = padic_radix_set_fmpz(t, cz, ctx);
        GR_IGNORE(padic_radix_set(r1, r0, ctx));
        if (s1 == GR_SUCCESS)
            s1 = padic_radix_submul(r1, x, t, ctx);
        GR_IGNORE(padic_radix_set(r2, r0, ctx));
        s2 = padic_radix_submul_fmpz(r2, x, cz, ctx);
        CHECK("submul_fmpz", r1, r2);

        /* product with an inexact zero: 0 + O(p^N) has valuation >= N */
        if (x->u.size == 0 && x->N != PADIC_RADIX_EXACT && !(y->u.size == 0 && y->N == PADIC_RADIX_EXACT))
        {
            slong vy = (y->u.size == 0) ? y->N : y->v;
            s1 = GR_SUCCESS;
            s2 = padic_radix_mul(r2, x, y, ctx);
            if (s2 == GR_SUCCESS && (r2->u.size != 0 ||
                    r2->N > FLINT_MIN(vy + x->N, PADIC_RADIX_ERR_MAX)))
            {
                flint_printf("FAIL: product with an inexact zero\n");
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }
        }

        /* pow vs one rounding of the exact power */
        {
            ulong n = n_randint(state, 10);

            s2 = padic_radix_pow_ui(r2, x, n, ctx);

            if (n >= 2 && x->u.size != 0 && FLINT_ABS(x->v) < 1000
                && (x->N == PADIC_RADIX_EXACT || FLINT_ABS(x->N) < 100000))
            {
                gr_ctx_t E;
                padic_radix_t X, P;
                fmpz_t nz, pz, rz;

                gr_ctx_init_padic_radix(E, p, PADIC_RADIX_PREC_INF, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
                padic_radix_init(X, E);
                padic_radix_init(P, E);
                radix_integer_set(&X->u, &x->u, PADIC_RADIX_CTX_RADIX(E));
                X->v = x->v;
                X->N = PADIC_RADIX_EXACT;
                GR_IGNORE(padic_radix_one(P, E));
                for (i = 0; i < (slong) n; i++)
                    GR_IGNORE(padic_radix_mul(P, P, X, E));

                fmpz_init_set_ui(nz, n);
                fmpz_init_set_ui(pz, p);
                fmpz_init(rz);
                P->N = (x->N == PADIC_RADIX_EXACT) ? PADIC_RADIX_EXACT
                     : (slong) (n - 1) * x->v + x->N + fmpz_remove(rz, nz, pz);
                s1 = _padic_radix_finalize(P, ctx);
                GR_IGNORE(padic_radix_set(r1, P, ctx));
                CHECK("pow_ui", r1, r2);

                fmpz_clear(nz);
                fmpz_clear(pz);
                fmpz_clear(rz);
                padic_radix_clear(X, E);
                padic_radix_clear(P, E);
                gr_ctx_clear(E);
            }

            fmpz_set_ui(cz, n);
            s1 = padic_radix_pow_fmpz(r1, x, cz, ctx);
            CHECK("pow_fmpz vs pow_ui", r1, r2);

            GR_IGNORE(padic_radix_set(r1, x, ctx));
            s1 = padic_radix_pow_ui(r1, r1, n, ctx);
            CHECK("pow_ui aliased", r1, r2);

            s1 = padic_radix_inv(r1, x, ctx);
            if (s1 == GR_SUCCESS)
                s1 = padic_radix_pow_ui(r1, r1, n, ctx);
            s2 = padic_radix_pow_si(r2, x, -(slong) n, ctx);
            if (n != 0)
                CHECK("pow_si vs inv + pow_ui", r1, r2);
        }

        /* vector-scalar operations vs elementwise loops */
        {
            slong len = n_randint(state, 10);
            padic_radix_struct * a, * v1, * v2;
            int op = n_randint(state, 6);

            a = gr_heap_init_vec(len, ctx);
            v1 = gr_heap_init_vec(len, ctx);
            v2 = gr_heap_init_vec(len, ctx);

            for (i = 0; i < len; i++)
            {
                _randelem(a + i, state, ctx);
                _randelem(v1 + i, state, ctx);
                GR_IGNORE(padic_radix_set(v2 + i, v1 + i, ctx));
            }

            s1 = s2 = GR_SUCCESS;
            fmpz_set_si(cz, c);

            for (i = 0; i < len; i++)
            {
                switch (op)
                {
                    case 0: s1 |= padic_radix_addmul(v1 + i, a + i, y, ctx); break;
                    case 1: s1 |= padic_radix_submul(v1 + i, a + i, y, ctx); break;
                    case 2: s1 |= padic_radix_mul(v1 + i, a + i, y, ctx); break;
                    case 3: s1 |= padic_radix_addmul_si(v1 + i, a + i, c, ctx); break;
                    case 4: s1 |= padic_radix_mul_si(v1 + i, a + i, c, ctx); break;
                    default: s1 |= padic_radix_addmul_fmpz(v1 + i, a + i, cz, ctx); break;
                }
            }

            switch (op)
            {
                case 0: s2 = _padic_radix_vec_addmul_scalar(v2, a, len, y, ctx); break;
                case 1: s2 = _padic_radix_vec_submul_scalar(v2, a, len, y, ctx); break;
                case 2: s2 = _padic_radix_vec_mul_scalar(v2, a, len, y, ctx); break;
                case 3: s2 = _padic_radix_vec_addmul_scalar_si(v2, a, len, c, ctx); break;
                case 4: s2 = _padic_radix_vec_mul_scalar_si(v2, a, len, c, ctx); break;
                default: s2 = _padic_radix_vec_addmul_scalar_fmpz(v2, a, len, cz, ctx); break;
            }

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
                for (i = 0; i < len; i++)
                    if (!_same(v1 + i, v2 + i, ctx))
                    {
                        flint_printf("FAIL: vector-scalar operation %d, index %wd\n", op, i);
                        gr_ctx_println(ctx);
                        flint_abort();
                    }

            gr_heap_clear_vec(a, len, ctx);
            gr_heap_clear_vec(v1, len, ctx);
            gr_heap_clear_vec(v2, len, ctx);
        }

        padic_radix_clear(x, ctx);
        padic_radix_clear(y, ctx);
        padic_radix_clear(r0, ctx);
        padic_radix_clear(r1, ctx);
        padic_radix_clear(r2, ctx);
        padic_radix_clear(t, ctx);
        fmpz_clear(cz);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
