/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "dfloat.h"

/* The vector functions and the scalar square roots with the output
   aliasing an input give bitwise the results of the out-of-place calls.
   (The SIMD blocks store their four results before redoing in scalar the
   lanes that need it -- for the ball square roots, every block mixing
   exact and inexact inputs -- and the scalar redo must not read the
   overwritten inputs; the plain square root kernels write the output
   before they have read all of the input.) */

#define ALIAS_UNARY(T, X, f) \
    do { \
        memcpy(r2, x, len * sizeof(T)); \
        (void) _ ## X ## _vec_ ## f((T *) r1, (const T *) x, len); \
        (void) _ ## X ## _vec_ ## f((T *) r2, (const T *) r2, len); \
        if (memcmp(r1, r2, len * sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X " vec_" #f " in place, len = %wd\n", len); \
    } while (0)

#define ALIAS_BINARY(T, X, f) \
    do { \
        (void) _ ## X ## _vec_ ## f((T *) r1, (const T *) x, (const T *) y, len); \
        memcpy(r2, x, len * sizeof(T)); \
        (void) _ ## X ## _vec_ ## f((T *) r2, (const T *) r2, (const T *) y, len); \
        if (memcmp(r1, r2, len * sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X " vec_" #f " in place (x), len = %wd\n", len); \
        memcpy(r2, y, len * sizeof(T)); \
        (void) _ ## X ## _vec_ ## f((T *) r2, (const T *) x, (const T *) r2, len); \
        if (memcmp(r1, r2, len * sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X " vec_" #f " in place (y), len = %wd\n", len); \
    } while (0)

#define ALIAS_SIN_COS(T, X) \
    do { \
        _ ## X ## _vec_sin_cos((T *) r1, (T *) s1, (const T *) x, len); \
        memcpy(r2, x, len * sizeof(T)); \
        _ ## X ## _vec_sin_cos((T *) r2, (T *) s2, (const T *) r2, len); \
        if (memcmp(r1, r2, len * sizeof(T)) != 0 || memcmp(s1, s2, len * sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X " vec_sin_cos in place (sin), len = %wd\n", len); \
        memcpy(s2, x, len * sizeof(T)); \
        _ ## X ## _vec_sin_cos((T *) r2, (T *) s2, (const T *) s2, len); \
        if (memcmp(r1, r2, len * sizeof(T)) != 0 || memcmp(s1, s2, len * sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X " vec_sin_cos in place (cos), len = %wd\n", len); \
    } while (0)

#define ALIAS_SCALAR_ROOT(T, X, f) \
    do { \
        T u, v; \
        memcpy(&u, x, sizeof(T)); \
        (void) X ## _ ## f(&v, &u); \
        (void) X ## _ ## f(&u, &u); \
        if (memcmp(&u, &v, sizeof(T)) != 0) \
            TEST_FUNCTION_FAIL(#X "_" #f " in place\n"); \
    } while (0)

#define ALIAS_ALL(T, X) \
    do { \
        switch (which) \
        { \
            case 0: ALIAS_UNARY(T, X, exp); break; \
            case 1: ALIAS_UNARY(T, X, sin); break; \
            case 2: ALIAS_UNARY(T, X, cos); break; \
            case 3: ALIAS_UNARY(T, X, sqrt); ALIAS_SCALAR_ROOT(T, X, sqrt); break; \
            case 4: ALIAS_UNARY(T, X, rsqrt); ALIAS_SCALAR_ROOT(T, X, rsqrt); break; \
            case 5: ALIAS_UNARY(T, X, log); break; \
            case 6: ALIAS_UNARY(T, X, atan); break; \
            case 7: ALIAS_SIN_COS(T, X); break; \
            case 8: ALIAS_BINARY(T, X, div); break; \
            case 9: ALIAS_BINARY(T, X, mul); break; \
            default: ALIAS_BINARY(T, X, add); break; \
        } \
    } while (0)

TEST_FUNCTION_START(vec_aliasing, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int n = 1 + n_randint(state, DFLOAT_MAX_N), ball = n_randint(state, 2);
        int which = n_randint(state, 11), mode = n_randint(state, 3);
        slong len = 1 + n_randint(state, 24), i, sz = 8 * (n + ball);
        double * x, * y, * r1, * r2, * s1, * s2;

        /* the partial functions (square roots, log, division) leave the
           elements where they fail unspecified: random special inputs
           only for the total ones */
        if (which == 3 || which == 4 || which == 5 || which == 8)
            mode = 0;

        x = flint_malloc(len * sz); y = flint_malloc(len * sz);
        r1 = flint_malloc(len * sz); r2 = flint_malloc(len * sz);
        s1 = flint_malloc(len * sz); s2 = flint_malloc(len * sz);

        /* positive values of moderate size (the fast paths), with the
           exact and inexact elements mixed within the SIMD blocks, or
           random special values */
        for (i = 0; i < len; i++)
        {
            double * xi = (double *) ((char *) x + i * sz), * yi = (double *) ((char *) y + i * sz);
            int k;
            if (mode == 0)
            {
                xi[0] = 0.25 + n_randint(state, 1000) / 128.0;
                yi[0] = 0.5 + n_randint(state, 1000) / 64.0;
                for (k = 1; k < n; k++)
                    xi[k] = yi[k] = 0.0;
                if (n > 1 && n_randint(state, 2))
                    xi[1] = ldexp(xi[0], -60) * (n_randint(state, 2) ? 1 : -1);
                if (ball)
                {
                    xi[n] = n_randint(state, 2) ? 0.0 : ldexp(xi[0], -40 - 50 * (n - 1));
                    yi[n] = n_randint(state, 2) ? 0.0 : ldexp(yi[0], -40 - 50 * (n - 1));
                }
            }
            else
            {
                switch (n + 4 * ball)
                {
                    case 1: d1_randtest((d1_ptr) xi, state); d1_randtest((d1_ptr) yi, state); break;
                    case 2: d2_randtest((d2_ptr) xi, state); d2_randtest((d2_ptr) yi, state); break;
                    case 3: d3_randtest((d3_ptr) xi, state); d3_randtest((d3_ptr) yi, state); break;
                    case 4: d4_randtest((d4_ptr) xi, state); d4_randtest((d4_ptr) yi, state); break;
                    case 5: d1b_randtest((d1b_ptr) xi, state); d1b_randtest((d1b_ptr) yi, state); break;
                    case 6: d2b_randtest((d2b_ptr) xi, state); d2b_randtest((d2b_ptr) yi, state); break;
                    case 7: d3b_randtest((d3b_ptr) xi, state); d3b_randtest((d3b_ptr) yi, state); break;
                    default: d4b_randtest((d4b_ptr) xi, state); d4b_randtest((d4b_ptr) yi, state); break;
                }
            }
        }

        switch (n + 4 * ball)
        {
            case 1: ALIAS_ALL(d1_struct, d1); break;
            case 2: ALIAS_ALL(d2_struct, d2); break;
            case 3: ALIAS_ALL(d3_struct, d3); break;
            case 4: ALIAS_ALL(d4_struct, d4); break;
            case 5: ALIAS_ALL(d1b_struct, d1b); break;
            case 6: ALIAS_ALL(d2b_struct, d2b); break;
            case 7: ALIAS_ALL(d3b_struct, d3b); break;
            default: ALIAS_ALL(d4b_struct, d4b); break;
        }

        flint_free(x); flint_free(y); flint_free(r1); flint_free(r2);
        flint_free(s1); flint_free(s2);
    }

    TEST_FUNCTION_END(state);
}
