/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_IMPL_H
#define FLINT_IMPL_H

/*
 * elementary functions
 */
/* NOTE: Casts to ulong. This is because -WORD_MIN is not defined as a slong. */
#define UI_ABS(x) ((slong)(x) < 0 ? -((ulong) (x)) : (ulong) (x))

#define r_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) >> (shift)))

#define l_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) << (shift)))



/*
 * swaps
 */
#define MP_PTR_SWAP(x, y)   \
    do {                    \
        mp_limb_t * __txxx; \
        __txxx = x;         \
        x = y;              \
        y = __txxx;         \
    } while (0)

#define SLONG_SWAP(A, B)    \
    do {                    \
        slong __t_m_p_ = A; \
        A = B;              \
        B = __t_m_p_;       \
    } while (0)

#define ULONG_SWAP(A, B)    \
    do {                    \
        ulong __t_m_p_ = A; \
        A = B;              \
        B = __t_m_p_;       \
    } while (0)

#define MP_LIMB_SWAP(A, B)      \
    do {                        \
        mp_limb_t __t_m_p_ = A; \
        A = B;                  \
        B = __t_m_p_;           \
    } while (0)

#define DOUBLE_SWAP(A, B)    \
    do {                     \
        double __t_m_p_ = A; \
        A = B;               \
        B = __t_m_p_;        \
    } while (0)



/*
 * mpn macros
 */
#define flint_mpn_zero(xxx, nnn)                \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < (nnn); ixxx++)    \
            (xxx)[ixxx] = UWORD(0);             \
    } while (0)

#define flint_mpn_copyi(xxx, yyy, nnn)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < (nnn); ixxx++)    \
            (xxx)[ixxx] = (yyy)[ixxx];          \
    } while (0)

#define flint_mpn_copyd(xxx, yyy, nnn)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = nnn - 1; ixxx >= 0; ixxx--) \
            (xxx)[ixxx] = (yyy)[ixxx];          \
    } while (0)

#define flint_mpn_store(xxx, nnn, yyy)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < nnn; ixxx++)      \
        (xxx)[ixxx] = yyy;                      \
    } while (0)



/*
 * Newton iteration macros
 */
#define FLINT_NEWTON_INIT(from, to)                     \
    {                                                   \
        slong __steps[FLINT_BITS], __i, __from, __to;   \
        __steps[__i = 0] = __to = (to);                 \
        __from = (from);                                \
        while (__to > __from)                           \
            __steps[++__i] = (__to = (__to + 1) / 2);   \

#define FLINT_NEWTON_BASECASE(bc_to)                    \
        { slong bc_to = __to;

#define FLINT_NEWTON_END_BASECASE                       \
        }

#define FLINT_NEWTON_LOOP(step_from, step_to)           \
        {                                               \
            for (__i--; __i >= 0; __i--)                \
            {                                           \
                slong step_from = __steps[__i+1];       \
                slong step_to = __steps[__i];           \

#define FLINT_NEWTON_END_LOOP                           \
            }}

#define FLINT_NEWTON_END                                \
    }



/*
 * garbage collection and test initializer and cleanup
 */
#if FLINT_USES_GC
#define FLINT_GC_INIT() GC_init()
#else
#define FLINT_GC_INIT()
#endif

#define FLINT_TEST_INIT(xxx)    \
   flint_rand_t xxx;            \
   FLINT_GC_INIT();             \
   flint_randinit(xxx)

#define FLINT_TEST_CLEANUP(xxx) \
   flint_randclear(xxx);        \
   flint_cleanup_master();


#endif
