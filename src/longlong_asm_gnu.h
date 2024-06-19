/*
    Copyright 1991-2018, 2021, 2022 Free Software Foundation, Inc.

    Copyright 2024 Albin Ahlbäck

    This file is part of FLINT.

    Contains code from GNU MP Library.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_ASM_GNU_H
#define LONGLONG_ASM_GNU_H

/* NOTE: Is missing from the immature CMake building script, so just try
 * something here. */
#ifndef LSYM_PREFIX
# define LSYM_PREFIX ".L"
#endif

/* ASM_L gives a local label for a gcc asm block, for use when temporary
   local labels like "1:" might not be available, which is the case for
   instance on the x86s (the SCO assembler doesn't support them).

   The label generated is made unique by including "%=" which is a unique
   number for each insn.  This ensures the same name can be used in multiple
   asm blocks, perhaps via a macro.  Since jumps between asm blocks are not
   allowed there's no need for a label to be usable outside a single
   block.  */

#define ASM_L(name)  LSYM_PREFIX "asm_%=_" #name

#if !FLINT_WANT_ASSERT && defined(__amd64__)
# define MPN_IORD_U(ptr, incr, aors) \
  do { \
    nn_ptr  __ptr_dummy; \
    if (__builtin_constant_p(incr) && (incr) == 0) \
    { \
    } \
    else if (__builtin_constant_p(incr) && (incr) == 1) \
    { \
      __asm__ volatile \
        ("\n" ASM_L(top) ":\n" \
         "\t" aors "\t$1, (%0)\n" \
         "\tlea\t%c2(%0), %0\n" \
         "\tjc\t" ASM_L(top) \
         : "=r" (__ptr_dummy) \
         : "0"  (ptr), "n" (sizeof(ulong)) \
         : "memory"); \
    } \
    else \
    { \
      __asm__ volatile \
        ( aors "\t%2, (%0)\n" \
         "\tjnc\t" ASM_L(done) "\n" \
         ASM_L(top) ":\n" \
         "\t" aors "\t$1, %c3(%0)\n" \
         "\tlea\t%c3(%0), %0\n" \
         "\tjc\t" ASM_L(top) "\n" \
         ASM_L(done) ":\n" \
         : "=r" (__ptr_dummy) \
         : "0"  (ptr), \
           "re" ((ulong) (incr)), "n" (sizeof(ulong)) \
         : "memory"); \
    } \
  } while (0)

# if FLINT_BITS == 32
#  define MPN_INCR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "addl")
#  define MPN_DECR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "subl")
# else
#  define MPN_INCR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "addq")
#  define MPN_DECR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "subq")
# endif
#elif !FLINT_WANT_ASSERT && (FLINT_BITS == 64 && defined(__aarch64__))
# define MPN_IORD_U(ptr, incr, aors, cond, anticond) \
  do { \
    nn_ptr  __ptr_dummy; \
    ulong  __reg_dummy; \
    if (__builtin_constant_p (incr) && (incr) == 0) \
    { \
    } \
    else if (__builtin_constant_p (incr) && (incr) == 1) \
    { \
      __asm__ volatile \
        ("\n" ASM_L(top) ":\n" \
         "\tldr\t%1, [%0]\n" \
         "\t" aors "\t%1, %1, #1\n" \
         "\tstr\t%1, [%0], #8\n" \
         "\tb." cond "\t" ASM_L(top) "\n" \
         : "=r" (__ptr_dummy), "=&r" (__reg_dummy) \
         : "0"  (ptr) \
         : "memory"); \
    } \
    else \
    { \
      __asm__ volatile \
        ("\tldr\t%1, [%0]\n" \
         "\t" aors "\t%1, %1, %3\n" \
         "\tstr\t%1, [%0], #8\n" \
         "\tb." anticond "\t" ASM_L(done) "\n" \
         "\n" ASM_L(top) ":\n" \
         "\tldr\t%1, [%0]\n" \
         "\t" aors "\t%1, %1, #1\n" \
         "\tstr\t%1, [%0], #8\n" \
         "\tb." cond "\t" ASM_L(top) "\n" \
         ASM_L(done) ":\n" \
         : "=r" (__ptr_dummy), "=&r" (__reg_dummy) \
         : "0"  (ptr), "rI" ((ulong) (incr)) \
         : "memory"); \
    } \
  } while (0)

# define MPN_INCR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "adds", "cs", "cc")
# define MPN_DECR_U(ptr, size, incr)  MPN_IORD_U(ptr, incr, "subs", "cc", "cs")
#endif

#endif
