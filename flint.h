/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_H
#define FLINT_H

#undef ulong
#define ulong ulongxx /* ensure vendor doesn't typedef ulong */
#if !defined(_MSC_VER)
#include <sys/param.h> /* for BSD define */
#endif
#include <gmp.h>
#include <mpfr.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h> /* for alloca on FreeBSD */
#if (!defined(BSD) && !defined(__MINGW64__) && !defined(__MINGW32__) && !defined(_MSC_VER)) || defined(__GNU__)
/* MinGW and FreeBSD have alloca, but not alloca.h */
#include <alloca.h>
#endif
#if defined(__MINGW32__)
#include <malloc.h> /* for alloca on MinGW */
#endif
#include "limits.h"
#include "longlong.h"
#include "flint-config.h"
#undef ulong

#ifdef FLINT_INLINES_C
#define FLINT_INLINE FLINT_DLL
#else
#define FLINT_INLINE static __inline__
#endif

#if FLINT_USES_GC
#include "gc.h"
#endif

#if FLINT_WANT_ASSERT
#include <assert.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

/* flint version number */

#define __FLINT_VERSION 2
#define __FLINT_VERSION_MINOR 9
#define __FLINT_VERSION_PATCHLEVEL 0
#define FLINT_VERSION "2.9.0"
#define __FLINT_RELEASE (__FLINT_VERSION * 10000 + \
                         __FLINT_VERSION_MINOR * 100 + \
                         __FLINT_VERSION_PATCHLEVEL)

/*
   Check mpir and mpfr version numbers
*/
#if __GNU_MP_VERSION < 5
#error GMP 5.0.0 or MPIR 2.6.0 or later are required
#endif

#if MPFR_VERSION_MAJOR < 3
#error MPFR 3.0.0 or later is required
#endif

/*
   We define alternative key words for "asm" and "inline", allowing
   the code to be compiled with the "-ansi" flag under GCC
 */
#ifndef __GNUC__
    #define __asm__     asm
    #define __inline__  inline
#endif

extern char flint_version[];

#define ulong mp_limb_t
#define slong mp_limb_signed_t

FLINT_DLL void * flint_malloc(size_t size);
FLINT_DLL void * flint_realloc(void * ptr, size_t size);
FLINT_DLL void * flint_calloc(size_t num, size_t size);
FLINT_DLL void flint_free(void * ptr);

typedef void (*flint_cleanup_function_t)(void);
FLINT_DLL void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function);
FLINT_DLL void flint_cleanup(void);
FLINT_DLL void flint_cleanup_master(void);

FLINT_DLL void __flint_set_memory_functions(void *(*alloc_func) (size_t),
     void *(*calloc_func) (size_t, size_t), void *(*realloc_func) (void *, size_t),
                                                              void (*free_func) (void *));

FLINT_DLL void __flint_get_memory_functions(void *(**alloc_func) (size_t),
     void *(**calloc_func) (size_t, size_t), void *(**realloc_func) (void *, size_t),
                                                              void (**free_func) (void *));

#ifdef __GNUC__
#define FLINT_NORETURN __attribute__ ((noreturn))
#else
#define FLINT_NORETURN
#endif

FLINT_DLL FLINT_NORETURN void flint_abort(void);
FLINT_DLL void flint_set_abort(FLINT_NORETURN void (*func)(void));
  /* flint_abort is calling abort by default
   * if flint_set_abort is used, then instead of abort this function
   * is called. EXPERIMENTALLY use at your own risk!
   * May disappear in future versions.
   */


#if defined(_WIN64) || defined(__mips64)
#if defined(__MINGW64__)
#define WORD_FMT "%I64"
#define WORD_WIDTH_FMT "%*I64"
#else
#define WORD_FMT "%ll"
#define WORD_WIDTH_FMT "%*ll"
#endif
#define WORD(xx) (xx##LL)
#define UWORD(xx) (xx##ULL)
#ifndef FLINT_NO_WORDMAC
#define UWORD_MAX ULLONG_MAX
#define UWORD_MIN ULLONG_MIN
#define WORD_MAX LLONG_MAX
#define WORD_MIN LLONG_MIN
#endif
#else
#define WORD_FMT "%l"
#define WORD_WIDTH_FMT "%*l"
#define WORD(xx) (xx##L)
#define UWORD(xx) (xx##UL)
#ifndef FLINT_NO_WORDMAC
#define UWORD_MAX ULONG_MAX
#define UWORD_MIN ULONG_MIN
#define WORD_MAX LONG_MAX
#define WORD_MIN LONG_MIN
#endif
#endif

#if GMP_LIMB_BITS == 64
    #define FLINT_BITS 64
    #define FLINT_D_BITS 53
    #define FLINT64 1
#else
    #define FLINT_BITS 32
    #define FLINT_D_BITS 31
#endif

#define flint_bitcnt_t ulong

#if FLINT_USES_TLS
#if defined(__GNUC__) && __STDC_VERSION__ >= 201112L && __GNUC__ == 4 && __GNUC_MINOR__ < 9
/* GCC 4.7, 4.8 with -std=gnu11 purport to support C11 via __STDC_VERSION__ but lack _Thread_local */
#define FLINT_TLS_PREFIX __thread
#elif __STDC_VERSION__ >= 201112L
#define FLINT_TLS_PREFIX _Thread_local
#elif defined(_MSC_VER)
#define FLINT_TLS_PREFIX __declspec(thread)
#elif defined(__GNUC__)
#define FLINT_TLS_PREFIX __thread
#else
#error "thread local prefix defined in C11 or later"
#endif
#else
#define FLINT_TLS_PREFIX
#endif

FLINT_DLL int flint_get_num_threads(void);
FLINT_DLL void flint_set_num_threads(int num_threads);
FLINT_DLL void _flint_set_num_workers(int num_workers);
FLINT_DLL int flint_set_num_workers(int num_workers);
FLINT_DLL void flint_reset_num_workers(int max_workers);
FLINT_DLL int flint_set_thread_affinity(int * cpus, slong length);
FLINT_DLL int flint_restore_thread_affinity();

int flint_test_multiplier(void);

typedef struct
{
    gmp_randstate_t gmp_state;
    int gmp_init;
    mp_limb_t __randval;
    mp_limb_t __randval2;
} flint_rand_s;

typedef flint_rand_s flint_rand_t[1];

FLINT_INLINE
void flint_randinit(flint_rand_t state)
{
   state->gmp_init = 0;
#if FLINT64
    state->__randval = UWORD(13845646450878251009);
    state->__randval2 = UWORD(13142370077570254774);
#else
    state->__randval = UWORD(4187301858);
    state->__randval2 = UWORD(3721271368);
#endif
}

FLINT_INLINE
void flint_randseed(flint_rand_t state, ulong seed1, ulong seed2)
{
   state->__randval = seed1;
   state->__randval2 = seed2;
}

FLINT_INLINE
void flint_get_randseed(ulong * seed1, ulong * seed2, flint_rand_t state)
{
   *seed1 = state->__randval;
   *seed2 = state->__randval2;
}


FLINT_INLINE
void _flint_rand_init_gmp(flint_rand_t state)
{
    if (!state->gmp_init)
    {
        gmp_randinit_default(state->gmp_state);
        state->gmp_init = 1;
    }
}

FLINT_INLINE
void flint_randclear(flint_rand_t state)
{
    if (state->gmp_init)
        gmp_randclear(state->gmp_state);
}

FLINT_INLINE
flint_rand_s * flint_rand_alloc(void)
{
    return (flint_rand_s *) flint_malloc(sizeof(flint_rand_s));
}

FLINT_INLINE
void flint_rand_free(flint_rand_s * state)
{
    flint_free(state);
}

#if FLINT_USES_GC
#define FLINT_GC_INIT() GC_init()
#else
#define FLINT_GC_INIT()
#endif

#define FLINT_TEST_INIT(xxx) \
   flint_rand_t xxx; \
   FLINT_GC_INIT(); \
   flint_randinit(xxx)

#define FLINT_TEST_CLEANUP(xxx) \
   flint_randclear(xxx); \
   flint_cleanup_master();

/*
  We define this here as there is no mpfr.h
 */
typedef __mpfr_struct flint_mpfr;

#if FLINT_WANT_ASSERT
#define FLINT_ASSERT(param) assert(param)
#else
#define FLINT_ASSERT(param)
#endif

#if defined(__GNUC__)
#define FLINT_UNUSED(x) UNUSED_ ## x __attribute__((unused))
#define FLINT_SET_BUT_UNUSED(x) x __attribute__((unused))
#if __GNUC__ >= 4
#define FLINT_WARN_UNUSED __attribute__((warn_unused_result))
#else
#define FLINT_WARN_UNUSED
#endif
#else
#define __attribute__(x)
#define FLINT_UNUSED(x) x
#define FLINT_SET_BUT_UNUSED(x) x
#define FLINT_WARN_UNUSED
#endif

#define FLINT_MAX(x, y) ((x) > (y) ? (x) : (y))
#define FLINT_MIN(x, y) ((x) > (y) ? (y) : (x))
#define FLINT_ABS(x) ((slong)(x) < 0 ? (-(x)) : (x))
#define FLINT_SIGN_EXT(x) (-(ulong)((slong)(x) < 0))
#define FLINT_SGN(x) ((0 < (slong)(x)) - ((slong)(x) < 0))

#define MP_PTR_SWAP(x, y) \
    do { \
        mp_limb_t * __txxx; \
        __txxx = x; \
        x = y; \
        y = __txxx; \
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

#define r_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) >> (shift)))

#define l_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) << (shift)))

#ifdef NEED_CLZ_TAB
FLINT_DLL extern const unsigned char __flint_clz_tab[128];
#endif

/* Beware when using the unsigned return value in signed arithmetic */
static __inline__
mp_limb_t FLINT_BIT_COUNT(mp_limb_t x)
{
   mp_limb_t zeros = FLINT_BITS;
   if (x) count_leading_zeros(zeros, x);
   return FLINT_BITS - zeros;
}

#define FLINT_FLOG2(k)  (FLINT_BIT_COUNT(k) - 1)

#define FLINT_CLOG2(k)  FLINT_BIT_COUNT((k) - 1)

#define flint_mpn_zero(xxx, nnn) \
    do \
    { \
        slong ixxx; \
        for (ixxx = 0; ixxx < (nnn); ixxx++) \
            (xxx)[ixxx] = UWORD(0); \
    } while (0)

#define flint_mpn_copyi(xxx, yyy, nnn) \
   do { \
      slong ixxx; \
      for (ixxx = 0; ixxx < (nnn); ixxx++) \
         (xxx)[ixxx] = (yyy)[ixxx]; \
   } while (0)

#define flint_mpn_copyd(xxx, yyy, nnn) \
   do { \
      slong ixxx; \
      for (ixxx = nnn - 1; ixxx >= 0; ixxx--) \
         (xxx)[ixxx] = (yyy)[ixxx]; \
   } while (0)

#define flint_mpn_store(xxx, nnn, yyy) \
   do \
   { \
      slong ixxx; \
      for (ixxx = 0; ixxx < nnn; ixxx++) \
         (xxx)[ixxx] = yyy; \
   } while (0)

/* common usage of flint_malloc */
#define FLINT_ARRAY_ALLOC(n, T) (T *) flint_malloc((n)*sizeof(T))
#define FLINT_ARRAY_REALLOC(p, n, T) (T *) flint_realloc(p, (n)*sizeof(T))

/* temporary allocation */
#define TMP_INIT \
   typedef struct __tmp_struct { \
      void * block; \
      struct __tmp_struct * next; \
   } __tmp_t; \
   __tmp_t * __tmp_root; \
   __tmp_t * __tpx

#define TMP_START \
   __tmp_root = NULL

#if FLINT_WANT_ASSERT
#define TMP_ALLOC(size) \
   (__tpx = (__tmp_t *) alloca(sizeof(__tmp_t)), \
       __tpx->next = __tmp_root, \
       __tmp_root = __tpx, \
       __tpx->block = flint_malloc(size))
#else
#define TMP_ALLOC(size) \
   (((size) > 8192) ? \
      (__tpx = (__tmp_t *) alloca(sizeof(__tmp_t)), \
       __tpx->next = __tmp_root, \
       __tmp_root = __tpx, \
       __tpx->block = flint_malloc(size)) : \
      alloca(size))
#endif

#define TMP_ARRAY_ALLOC(n, T) (T *) TMP_ALLOC((n)*sizeof(T))

#define TMP_END \
   while (__tmp_root) { \
      flint_free(__tmp_root->block); \
      __tmp_root = __tmp_root->next; \
   }

/* compatibility between gmp and mpir */
#ifndef mpn_com_n
#define mpn_com_n mpn_com
#endif

#ifndef mpn_neg_n
#define mpn_neg_n mpn_neg
#endif

#ifndef mpn_tdiv_q
/* substitute for mpir's mpn_tdiv_q */
static __inline__ void
mpn_tdiv_q(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn)
{
    mp_ptr _scratch;
    TMP_INIT;
    TMP_START;
    _scratch = (mp_ptr) TMP_ALLOC(dn * sizeof(mp_limb_t));
    mpn_tdiv_qr(qp, _scratch, 0, np, nn, dp, dn);
    TMP_END;
}
#endif

/* Newton iteration macros */
#define FLINT_NEWTON_INIT(from, to) \
    { \
        slong __steps[FLINT_BITS], __i, __from, __to; \
        __steps[__i = 0] = __to = (to); \
        __from = (from); \
        while (__to > __from) \
            __steps[++__i] = (__to = (__to + 1) / 2); \

#define FLINT_NEWTON_BASECASE(bc_to) { slong bc_to = __to;

#define FLINT_NEWTON_END_BASECASE }

#define FLINT_NEWTON_LOOP(step_from, step_to) \
        { \
            for (__i--; __i >= 0; __i--) \
            { \
                slong step_from = __steps[__i+1]; \
                slong step_to = __steps[__i]; \

#define FLINT_NEWTON_END_LOOP }}

#define FLINT_NEWTON_END }

FLINT_DLL int parse_fmt(int * floating, const char * fmt);

FLINT_DLL int flint_printf(const char * str, ...); /* flint version of printf */
FLINT_DLL int flint_vprintf(const char * str, va_list ap); /* va_list version of flint_printf */
FLINT_DLL int flint_fprintf(FILE * f, const char * str, ...); /* flint version of fprintf */
FLINT_DLL int flint_sprintf(char * s, const char * str, ...); /* flint version of sprintf */

FLINT_DLL int flint_scanf(const char * str, ...); /* flint version of scanf */
FLINT_DLL int flint_fscanf(FILE * f, const char * str, ...); /* flint version of fscanf */
FLINT_DLL int flint_sscanf(const char * s, const char * str, ...); /* flint version of sscanf */

FLINT_INLINE slong flint_mul_sizes(slong x, slong y)
{
    ulong hi, lo;

    umul_ppmm(hi, lo, (ulong) x, (ulong) y);
    if (hi != 0 || lo > WORD_MAX)
    {
        flint_printf("Exception (flint). Overflow creating size %wd x %wd object.\n", x, y);
        flint_abort();
    }
    return lo;
}

#include "gmpcompat.h"
#include "exception.h"

#ifdef __cplusplus
}
#endif

#endif
