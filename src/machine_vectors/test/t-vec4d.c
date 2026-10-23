/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "machine_vectors.h"

/* (32-byte generic vectors on a target without AVX, e.g. 32-bit x86,
   are passed differently from AVX vectors; nothing here crosses a
   library boundary) */
#if defined(__clang__)
# if defined(__has_warning)
#  if __has_warning("-Wpsabi")
#   pragma clang diagnostic ignored "-Wpsabi"
#  endif
# endif
#elif defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wpsabi"
#endif

/* The vec4d masks, bitwise and integer-lane operations, fused
   multiply-adds, permutes and gathers (used by dfloat), lane by lane
   against scalar references, with the backend that machine_vectors.h
   selects for this build. */

static uint64_t
_t_bits(double x)
{
    uint64_t u;
    memcpy(&u, &x, sizeof(u));
    return u;
}

static double
_t_from_bits(uint64_t u)
{
    double x;
    memcpy(&x, &u, sizeof(x));
    return x;
}

static double
_t_rand_d(flint_rand_t state)
{
    switch (n_randint(state, 12))
    {
        case 0: return 0.0;
        case 1: return -0.0;
        case 2: return INFINITY;
        case 3: return -INFINITY;
        case 4: return NAN;
        case 5: return (double) ((slong) n_randint(state, 21) - 10);
        case 6: return ((double) ((slong) n_randint(state, 21) - 10)) + 0.5;
        case 7: return _t_from_bits(n_randlimb(state) | ((uint64_t) n_randlimb(state) << 32));
        default:
            return ldexp((double) n_randlimb(state) / 4294967296.0 - 0.5,
                (int) n_randint(state, 200) - 100);
    }
}

static int
_t_same(double a, double b)
{
    return _t_bits(a) == _t_bits(b);
}

static void
_t_get(double * r, vec4d a)
{
    vec4d_store_unaligned(r, a);
}

#define T_CHECK(cond, what) \
    do { \
        if (!(cond)) \
        { \
            flint_printf("FAIL: %s, lane %d\n", what, i); \
            flint_printf("x = %.17g %.17g %.17g %.17g\n", x[0], x[1], x[2], x[3]); \
            flint_printf("y = %.17g %.17g %.17g %.17g\n", y[0], y[1], y[2], y[3]); \
            flint_abort(); \
        } \
    } while (0)

/* whether the C library's fma is correctly rounded (it is emulated,
   sometimes incorrectly, on some platforms without hardware fma) */
static int
_t_fma_ok(void)
{
    volatile double a = 1.0 + ldexp(1.0, -30), c = -(1.0 + ldexp(1.0, -29));
    return fma(a, a, c) == ldexp(1.0, -60);
}

TEST_FUNCTION_START(machine_vectors_vec4d, state)
{
    slong iter;
    int fma_ok = _t_fma_ok();

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        double x[4], y[4], z[4], r[4], t[20];
        uint64_t u[4], w[4];
        vec4d a, b, c;
        int i, m, n;

        for (i = 0; i < 4; i++)
        {
            x[i] = _t_rand_d(state);
            y[i] = n_randint(state, 4) ? _t_rand_d(state) : x[i];
            z[i] = _t_rand_d(state);
        }
        a = vec4d_load_unaligned(x);
        b = vec4d_load_unaligned(y);
        c = vec4d_load_unaligned(z);

        /* broadcast keeps the sign of zero */
        _t_get(r, vec4d_set_d(x[0]));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_same(r[i], x[0]), "set_d");
        _t_get(r, vec4d_set_d4(x[3], x[2], x[1], x[0]));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_same(r[i], x[3 - i]), "set_d4");

#define T_CMP(name, pred) \
        _t_get(r, vec4d_cmp_##name(a, b)); \
        for (i = 0; i < 4; i++) \
            T_CHECK(_t_bits(r[i]) == ((pred) ? ~UINT64_C(0) : UINT64_C(0)), "cmp_" #name);
        T_CMP(eq, x[i] == y[i])
        T_CMP(ne, !(x[i] == y[i]))
        T_CMP(lt, x[i] < y[i])
        T_CMP(le, x[i] <= y[i])
        T_CMP(gt, x[i] > y[i])
        T_CMP(ge, x[i] >= y[i])
        T_CMP(nlt, !(x[i] < y[i]))
        T_CMP(nle, !(x[i] <= y[i]))
        T_CMP(ngt, !(x[i] > y[i]))
        T_CMP(nge, !(x[i] >= y[i]))
#undef T_CMP

        /* blendv on the sign bit, with a mask or any value */
        _t_get(r, vec4d_blendv(a, b, c));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_same(r[i], (_t_bits(z[i]) >> 63) ? y[i] : x[i]), "blendv");
        _t_get(r, vec4d_blendv(a, b, vec4d_cmp_lt(a, b)));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_same(r[i], (x[i] < y[i]) ? y[i] : x[i]), "blendv mask");

        m = vec4d_movemask(c);
        for (i = 0; i < 4; i++)
            T_CHECK(((m >> i) & 1) == (int) (_t_bits(z[i]) >> 63), "movemask");
        T_CHECK((m >> 4) == 0, "movemask");

#define T_BIT(name, op) \
        _t_get(r, vec4d_bit_##name(a, b)); \
        for (i = 0; i < 4; i++) \
            T_CHECK(_t_bits(r[i]) == (op), "bit_" #name);
        T_BIT(and, _t_bits(x[i]) & _t_bits(y[i]))
        T_BIT(or, _t_bits(x[i]) | _t_bits(y[i]))
        T_BIT(xor, _t_bits(x[i]) ^ _t_bits(y[i]))
        T_BIT(andnot, ~_t_bits(x[i]) & _t_bits(y[i]))
#undef T_BIT

        for (i = 0; i < 4; i++)
        {
            u[i] = _t_bits(x[i]);
            w[i] = _t_bits(y[i]);
        }
        n = (int) n_randint(state, 64);
        _t_get(r, vec4d_bits_add(a, b));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_bits(r[i]) == u[i] + w[i], "bits_add");
        _t_get(r, vec4d_bits_sub(a, b));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_bits(r[i]) == u[i] - w[i], "bits_sub");
        _t_get(r, vec4d_bits_shl(a, n));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_bits(r[i]) == (u[i] << n), "bits_shl");
        _t_get(r, vec4d_bits_shr(a, n));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_bits(r[i]) == (u[i] >> n), "bits_shr");
        _t_get(r, vec4d_bits_set_u64(u[1]));
        for (i = 0; i < 4; i++)
            T_CHECK(_t_bits(r[i]) == u[1], "bits_set_u64");

        /* correctly rounded operations (nan payloads are not compared) */
#define T_FP(what, expr, ref) \
        _t_get(r, expr); \
        for (i = 0; i < 4; i++) \
        { \
            double e = (ref); \
            T_CHECK(isnan(e) ? isnan(r[i]) : _t_same(r[i], e), what); \
        }
        T_FP("sqrt", vec4d_sqrt(a), sqrt(x[i]))
        /* (against the C library's fma, if that is correctly rounded) */
        if (fma_ok)
        {
            T_FP("fmadd_fused", vec4d_fmadd_fused(a, b, c), fma(x[i], y[i], z[i]))
            T_FP("fmsub_fused", vec4d_fmsub_fused(a, b, c), fma(x[i], y[i], -z[i]))
            T_FP("fnmadd_fused", vec4d_fnmadd_fused(a, b, c), fma(-x[i], y[i], z[i]))
        }
        T_FP("div", vec4d_div(a, b), x[i] / y[i])
        T_FP("abs", vec4d_abs(a), fabs(x[i]))
        T_FP("neg", vec4d_neg(a), -x[i])
        T_FP("round", vec4d_round(a), nearbyint(x[i]))
        T_FP("floor", vec4d_floor(a), floor(x[i]))
#undef T_FP

        /* min and max away from nans, where the backends differ */
        if (!isnan(x[0]) && !isnan(x[1]) && !isnan(x[2]) && !isnan(x[3]) &&
            !isnan(y[0]) && !isnan(y[1]) && !isnan(y[2]) && !isnan(y[3]))
        {
            _t_get(r, vec4d_min(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(r[i] == FLINT_MIN(x[i], y[i]), "min");
            _t_get(r, vec4d_max(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(r[i] == FLINT_MAX(x[i], y[i]), "max");
        }

        /* permutes */
#define T_PERM(i0, i1, i2, i3) \
        _t_get(r, vec4d_permute_##i0##_##i1##_##i2##_##i3(a)); \
        { \
            const int p[4] = {i0, i1, i2, i3}; \
            for (i = 0; i < 4; i++) \
                T_CHECK(_t_same(r[i], x[p[i]]), "permute"); \
        }
        T_PERM(0, 1, 1, 0)
        T_PERM(0, 2, 0, 2)
        T_PERM(0, 2, 1, 3)
        T_PERM(1, 3, 1, 3)
        T_PERM(3, 1, 2, 0)
        T_PERM(3, 2, 1, 0)
#undef T_PERM
        {
            double e0[4] = {x[0], x[1], y[0], y[1]};
            double e1[4] = {x[2], x[3], y[2], y[3]};
            double e2[4] = {x[0], y[0], x[2], y[2]};
            double e3[4] = {x[1], y[1], x[3], y[3]};
            _t_get(r, vec4d_permute2_0_2(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(_t_same(r[i], e0[i]), "permute2_0_2");
            _t_get(r, vec4d_permute2_1_3(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(_t_same(r[i], e1[i]), "permute2_1_3");
            _t_get(r, vec4d_unpacklo(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(_t_same(r[i], e2[i]), "unpacklo");
            _t_get(r, vec4d_unpackhi(a, b));
            for (i = 0; i < 4; i++)
                T_CHECK(_t_same(r[i], e3[i]), "unpackhi");
        }

        /* gathers from integer-valued lanes */
        {
            double idx[4];
            for (i = 0; i < 20; i++)
                t[i] = _t_rand_d(state);
            for (i = 0; i < 4; i++)
                idx[i] = (double) n_randint(state, 20);
            _t_get(r, vec4d_gather(t, vec4d_to_index(vec4d_load_unaligned(idx))));
            for (i = 0; i < 4; i++)
                T_CHECK(_t_same(r[i], t[(int) idx[i]]), "gather");
        }
    }

    TEST_FUNCTION_END(state);
}
#undef T_CHECK
