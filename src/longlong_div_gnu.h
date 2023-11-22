/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_DIV_H
#define LONGLONG_DIV_H

#if defined (__amd64__) || (GMP_LIMB_BITS == 32 && (defined (__i386__) || defined (__i486__)))

# if GMP_LIMB_BITS == 64 && defined (__amd64__)
#  define _FLINT_ASM_DIV "divq"
#  define _FLINT_ASM_IDIV "idivq"
# else
#  define _FLINT_ASM_DIV "divl"
#  define _FLINT_ASM_IDIV "idivl"
# endif

# define udiv_qrnnd(q, r, n1, n0, dx) \
  __asm__ volatile (_FLINT_ASM_DIV " %4" \
    : "=a" (q), "=d" (r) \
    : "0" ((ulong)(n0)), "1" ((ulong)(n1)), \
      "rm" ((ulong)(dx)))

# define sdiv_qrnnd(q, r, n1, n0, dx) \
  __asm__ volatile (_FLINT_ASM_IDIV " %4" \
    : "=a" (q), "=d" (r) \
    : "0" ((ulong)(n0)), "1" ((ulong)(n1)), \
      "rm" ((ulong)(dx)))
#endif

#endif
