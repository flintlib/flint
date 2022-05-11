/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MINI_H
#define NMOD_MINI_H

#ifdef NMOD_INLINES_C
#define NMOD_INLINE FLINT_DLL
#else
#define NMOD_INLINE static __inline__
#endif

#include "ulong_extras_mini.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define NMOD_RED(r, a, mod)             \
     do                                 \
     {                                  \
         NMOD_RED2(r, UWORD(0), a, mod);\
     } while (0)

#define NMOD_RED2(r, a_hi, a_lo, mod)               \
    do                                              \
    {                                               \
        ulong q0xx, q1xx, r1xx;                     \
        const ulong u1xx = ((a_hi)<<(mod).norm)     \
            + ((mod.norm == 0)                      \
                    ? UWORD(0)                      \
                    : (ulong) (a_lo) >> (FLINT_BITS - (mod).norm)); \
        const ulong u0xx = ((a_lo)<<(mod).norm);    \
        const ulong nxx = ((mod).n<<(mod).norm);    \
        umul_ppmm(q1xx, q0xx, (mod).ninv, u1xx);    \
        add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx); \
        r1xx = (u0xx - (q1xx + 1)*nxx);             \
        if (r1xx > q0xx)                            \
            r1xx += nxx;                            \
        if (r1xx < nxx)                             \
            r = (r1xx>>(mod).norm);                 \
        else                                        \
            r = ((r1xx - nxx)>>(mod).norm);         \
    } while (0)

#define NMOD_RED3(r, a_hi, a_me, a_lo, mod) \
     do                                     \
     {                                      \
         ulong v_hi;	                    \
         NMOD_RED2(v_hi, a_hi, a_me, mod);  \
         NMOD_RED2(r, v_hi, a_lo, mod);     \
     } while (0)

#define NMOD2_RED2(r, a_hi, a_lo, mod)  \
     do                                 \
     {                                  \
         ulong v_hi;	                \
         NMOD_RED(v_hi, a_hi, mod);     \
         NMOD_RED2(r, v_hi, a_lo, mod); \
     } while (0)

#ifdef __cplusplus
}
#endif

#endif
