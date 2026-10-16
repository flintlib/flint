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
#include "ulong_extras.h"
#include "machine_vectors.h"

#if defined(__clang__)
# if defined(__has_warning)
#  if __has_warning("-Wpsabi")
#   pragma clang diagnostic ignored "-Wpsabi"
#  endif
# endif
#elif defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wpsabi"
#endif

/* The vec4n operations used by the radix-2^56 shift-and-add loops of
   mp_real (loads and stores, add, sub, bitwise and/or, shifts by a
   common count, the lane shift), lane by lane against scalar
   references, with the backend that machine_vectors.h selects for
   this build.  Only meaningful with 64-bit lanes. */

#define T_CHECK(name, z, ref) \
    do { \
        ulong _o[4]; \
        int _i; \
        vec4n_store_unaligned(_o, z); \
        for (_i = 0; _i < 4; _i++) \
            if (_o[_i] != (ref)[_i]) \
                TEST_FUNCTION_FAIL("%s: lane %d: %wu != %wu\n", name, _i, \
                    _o[_i], (ref)[_i]); \
    } while (0)

TEST_FUNCTION_START(machine_vectors_vec4n, state)
{
#if FLINT_BITS == 64
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong a[4], b[4], al[8], ref[4];
        vec4n va, vb;
        ulong s = n_randint(state, FLINT_BITS);
        int i;

        /* aligned storage for the aligned load/store */
        ulong * ap = (ulong *) (((uintptr_t) al + 31) & ~(uintptr_t) 31);

        for (i = 0; i < 4; i++)
        {
            a[i] = n_randtest(state);
            b[i] = n_randtest(state);
        }
        va = vec4n_load_unaligned(a);
        vb = vec4n_load_unaligned(b);

        T_CHECK("load_unaligned", va, a);
        memcpy(ap, a, sizeof(a));
        T_CHECK("load_aligned", vec4n_load_aligned(ap), a);
        vec4n_store_aligned(ap, vb);
        T_CHECK("store_aligned", vec4n_load_unaligned(ap), b);

        for (i = 0; i < 4; i++) ref[i] = a[i] + b[i];
        T_CHECK("add", vec4n_add(va, vb), ref);
        for (i = 0; i < 4; i++) ref[i] = a[i] - b[i];
        T_CHECK("sub", vec4n_sub(va, vb), ref);
        for (i = 0; i < 4; i++) ref[i] = a[i] & b[i];
        T_CHECK("bit_and", vec4n_bit_and(va, vb), ref);
        for (i = 0; i < 4; i++) ref[i] = a[i] | b[i];
        T_CHECK("bit_or", vec4n_bit_or(va, vb), ref);
        for (i = 0; i < 4; i++) ref[i] = a[i] >> s;
        T_CHECK("bit_shift_right", vec4n_bit_shift_right(va, s), ref);
        for (i = 0; i < 4; i++) ref[i] = a[i] << s;
        T_CHECK("bit_shift_left", vec4n_bit_shift_left(va, s), ref);

        ref[0] = b[3]; ref[1] = a[0]; ref[2] = a[1]; ref[3] = a[2];
        T_CHECK("shift_lanes_up", vec4n_shift_lanes_up(va, vb), ref);

        for (i = 0; i < 4; i++) ref[i] = a[0];
        T_CHECK("set_n", vec4n_set_n(a[0]), ref);
        for (i = 0; i < 4; i++) ref[i] = 0;
        T_CHECK("zero", vec4n_zero(), ref);
    }
#endif

    TEST_FUNCTION_END(state);
}
