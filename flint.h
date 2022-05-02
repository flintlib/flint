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
   Check GMP/MPIR version numbers
*/
#if __GNU_MP_VERSION < 5
#error GMP 5.0.0 or MPIR 2.6.0 or later are required
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

#if GMP_LIMB_BITS == 64
    #define FLINT_BITS 64
    #define FLINT_D_BITS 53
    #define FLINT64 1
#else
    #define FLINT_BITS 32
    #define FLINT_D_BITS 31
#endif

/*##############################################################################
# external type definitions
##############################################################################*/

#define ulong mp_limb_t
#define slong mp_limb_signed_t

#define flint_bitcnt_t ulong

typedef struct
{
    gmp_randstate_t gmp_state;
    int gmp_init;
    mp_limb_t __randval;
    mp_limb_t __randval2;
} flint_rand_s;

typedef flint_rand_s flint_rand_t[1];

typedef struct
{
    double * entries;
    slong r;
    slong c;
    double ** rows;
} d_mat_struct;

typedef d_mat_struct d_mat_t[1];

typedef struct
{
   mp_limb_t n;
   mp_limb_t ninv;
   flint_bitcnt_t norm;
} nmod_t;

typedef struct
{
    mp_ptr coeffs;
    slong alloc;
    slong length;
    nmod_t mod;
} nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

typedef struct
{
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    nmod_t mod;
}
nmod_mat_struct;

typedef nmod_mat_struct nmod_mat_t[1];

typedef slong fmpz;
typedef fmpz fmpz_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
} fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

typedef struct
{
    fmpz * entries;
    slong r;
    slong c;
    fmpz ** rows;
} fmpz_mat_struct;

typedef fmpz_mat_struct fmpz_mat_t[1];

typedef struct fmpz_mod_ctx {
    fmpz_t n;
    void (* add_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* sub_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* mul_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    nmod_t mod;
    ulong n_limbs[3];
    ulong ninv_limbs[3];
} fmpz_mod_ctx_struct;

typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_poly_struct;

typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

typedef struct
{
    fmpz_mat_t mat;
    fmpz_t mod;
}
fmpz_mod_mat_struct;

typedef fmpz_mod_mat_struct fmpz_mod_mat_t[1];

typedef struct
{
    fmpz num;
    fmpz den;
}
fmpq;

typedef fmpq fmpq_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    fmpz_t den;
} fmpq_poly_struct;

typedef fmpq_poly_struct fmpq_poly_t[1];

typedef struct
{
    fmpq * entries;
    slong r;
    slong c;
    fmpq ** rows;
} fmpq_mat_struct;

typedef fmpq_mat_struct fmpq_mat_t[1];

typedef struct
{
    fmpz_poly_struct * num;
    fmpz_poly_struct * den;
} fmpz_poly_q_struct;

typedef fmpz_poly_q_struct fmpz_poly_q_t[1];

typedef struct
{
    fmpz u;
    slong v;
    slong N;
} padic_struct;

typedef padic_struct padic_t[1];

enum padic_print_mode
{
    PADIC_TERSE, 
    PADIC_SERIES, 
    PADIC_VAL_UNIT
};

typedef struct
{
    fmpz_t p;
    double pinv;
    fmpz * pow;
    slong min;
    slong max;
    enum padic_print_mode mode;
} padic_ctx_struct;

typedef padic_ctx_struct padic_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    slong val;
    slong N;
} padic_poly_struct;

typedef padic_poly_struct padic_poly_t[1];

typedef struct
{
    fmpz_mat_struct mat;
    slong val;
    slong N;
} padic_mat_struct;

typedef padic_mat_struct padic_mat_t[1];

typedef padic_poly_struct qadic_struct;
typedef padic_poly_t qadic_t;

typedef struct
{
    padic_ctx_struct pctx;
    fmpz * a;
    slong * j;
    slong len;
    char * var;
} qadic_ctx_struct;

typedef qadic_ctx_struct qadic_ctx_t[1];

typedef fmpz_poly_t fq_t;
typedef fmpz_poly_struct fq_struct;

typedef struct
{
    fmpz_mod_ctx_t ctxp;
    int sparse_modulus;
    int is_conway; /* whether field was initialized with the Flint Conway tables (assures primitivity) */
    fmpz * a;
    slong * j;
    slong len;
    fmpz_mod_poly_t modulus;
    fmpz_mod_poly_t inv;
    char *var;
} fq_ctx_struct;

typedef fq_ctx_struct fq_ctx_t[1];

typedef struct
{
    fq_struct * coeffs;
    slong alloc;
    slong length;
} fq_poly_struct;

typedef fq_poly_struct fq_poly_t[1];

typedef struct
{
    fq_struct * entries;
    slong r;
    slong c;
    fq_struct ** rows;
} fq_mat_struct;

typedef fq_mat_struct fq_mat_t[1];

typedef nmod_poly_t fq_nmod_t;
typedef nmod_poly_struct fq_nmod_struct;

typedef struct
{
    fmpz p;
    nmod_t mod;
    int sparse_modulus;
    int is_conway; /* whether field was generated using Flint Conway table (assures primitivity) */
    mp_limb_t * a;
    slong * j;
    slong len;
    nmod_poly_t modulus;
    nmod_poly_t inv;
    char * var;
} fq_nmod_ctx_struct;

typedef fq_nmod_ctx_struct fq_nmod_ctx_t[1];

typedef struct
{
    fq_nmod_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_poly_struct;

typedef fq_nmod_poly_struct fq_nmod_poly_t[1];

typedef struct
{
    fq_nmod_struct * entries;
    slong r;
    slong c;
    fq_nmod_struct ** rows;
} fq_nmod_mat_struct;

typedef fq_nmod_mat_struct fq_nmod_mat_t[1];

typedef struct
{
    mp_limb_t value;
} fq_zech_struct;

typedef fq_zech_struct fq_zech_t[1];

typedef struct
{
    mp_limb_t qm1;              /* q - 1 */
    mp_limb_t qm1o2;            /* (q - 1) / 2 or 1 when p == 2 */
    mp_limb_t qm1opm1;          /* (q - 1) / (p - 1) */
    mp_limb_t p;
    double ppre;
    mp_limb_t prime_root;       /* primitive root for prime subfield */
    mp_limb_t * zech_log_table;
    mp_limb_t * prime_field_table;
    mp_limb_t * eval_table;
    fq_nmod_ctx_struct * fq_nmod_ctx;
    int owns_fq_nmod_ctx;
    int is_conway; /* whether field was generated using Flint Conway tables (assures primitivity) */
} fq_zech_ctx_struct;

typedef fq_zech_ctx_struct fq_zech_ctx_t[1];

typedef struct
{
    fq_zech_struct * coeffs;
    slong alloc;
    slong length;
} fq_zech_poly_struct;

typedef fq_zech_poly_struct fq_zech_poly_t[1];

typedef struct
{
    fq_zech_struct * entries;
    slong r;
    slong c;
    fq_zech_struct ** rows;
} fq_zech_mat_struct;

typedef fq_zech_mat_struct fq_zech_mat_t[1];

typedef union fq_default_struct
{
    fq_t fq;
    fq_nmod_t fq_nmod;
    fq_zech_t fq_zech;
    ulong nmod;
    fmpz_t fmpz_mod;
} fq_default_struct;

typedef fq_default_struct fq_default_t[1];

typedef struct
{
    int type;
    union ctx
    {
        fq_ctx_t fq;
        fq_nmod_ctx_t fq_nmod;
        fq_zech_ctx_t fq_zech;
        struct {
            nmod_t mod;
            mp_limb_t a;    /* minpoly is x - a */
        } nmod;
        struct {
            fmpz_mod_ctx_t mod;
            fmpz_t a;       /* minpoly is x - a */
        } fmpz_mod;
    } ctx;
} fq_default_ctx_struct;

typedef fq_default_ctx_struct fq_default_ctx_t[1];

typedef union fq_default_poly_struct
{
    fq_poly_t fq;
    fq_nmod_poly_t fq_nmod;
    fq_zech_poly_t fq_zech;
    nmod_poly_t nmod;
    fmpz_mod_poly_t fmpz_mod;
}
fq_default_poly_struct;

typedef fq_default_poly_struct fq_default_poly_t[1];

typedef union fq_default_mat_struct
{
    fq_mat_t fq;
    fq_nmod_mat_t fq_nmod;
    fq_zech_mat_t fq_zech;
    nmod_mat_t nmod;
    fmpz_mod_mat_t fmpz_mod;
} fq_default_mat_struct;

typedef fq_default_mat_struct fq_default_mat_t[1];

typedef struct
{
    nmod_poly_struct * entries;
    slong r;
    slong c;
    nmod_poly_struct ** rows;
    mp_limb_t modulus;
} nmod_poly_mat_struct;

typedef nmod_poly_mat_struct nmod_poly_mat_t[1];

typedef struct
{
    fmpz_poly_struct * entries;
    slong r;
    slong c;
    fmpz_poly_struct ** rows;
} fmpz_poly_mat_struct;

typedef fmpz_poly_mat_struct fmpz_poly_mat_t[1];

typedef struct
{
    int sign;
    fmpz * p;
    ulong * exp;
    slong alloc;
    slong num;
} fmpz_factor_struct;

typedef fmpz_factor_struct fmpz_factor_t[1];

typedef struct
{
    nmod_poly_struct * p;
    slong * exp;
    slong num;
    slong alloc;
} nmod_poly_factor_struct;

typedef nmod_poly_factor_struct nmod_poly_factor_t[1];

typedef struct {
    fmpz c;
    fmpz_poly_struct *p;
    slong *exp;
    slong num;
    slong alloc;
} fmpz_poly_factor_struct;

typedef fmpz_poly_factor_struct fmpz_poly_factor_t[1];

typedef struct
{
    fmpz_mod_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fmpz_mod_poly_factor_struct;

typedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1];

typedef struct
{
    fq_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_poly_factor_struct;

typedef fq_poly_factor_struct fq_poly_factor_t[1];

typedef struct
{
    fq_nmod_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_nmod_poly_factor_struct;

typedef fq_nmod_poly_factor_struct fq_nmod_poly_factor_t[1];

typedef struct
{
    fq_zech_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_zech_poly_factor_struct;

typedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1];

typedef union fq_default_poly_factor_struct
{
    fq_poly_factor_t fq;
    fq_nmod_poly_factor_t fq_nmod;
    fq_zech_poly_factor_t fq_zech;
    nmod_poly_factor_t nmod;
    fmpz_mod_poly_factor_t fmpz_mod;
} fq_default_poly_factor_struct;

typedef fq_default_poly_factor_struct fq_default_poly_factor_t[1];

typedef enum {
    ORD_LEX,
    ORD_DEGLEX,
    ORD_DEGREVLEX
} ordering_t;

#define MPOLY_NUM_ORDERINGS 3

typedef struct
{
    slong nvars;    /* number of variables */
    slong nfields;  /* number of fields in exponent vector */
    ordering_t ord; /* monomial ordering */
    int deg;        /* is ord a degree ordering? */
    int rev;        /* is ord a reversed ordering? */
    slong lut_words_per_exp[FLINT_BITS];
    unsigned char lut_fix_bits[FLINT_BITS]; /* FLINT_BITS < 256 */
} mpoly_ctx_struct;

typedef mpoly_ctx_struct mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    nmod_t mod;
} nmod_mpoly_ctx_struct;

typedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1];

typedef struct
{
    mp_limb_t * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fmpz_mod_ctx_t ffinfo;
} fmpz_mod_mpoly_ctx_struct;

typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} fmpz_mod_mpoly_struct;

typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

typedef struct
{
    fmpz_mpoly_ctx_t zctx;
} fmpq_mpoly_ctx_struct;

typedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpoly_ctx_struct;

typedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1];

typedef struct {
    mp_limb_t * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} fq_nmod_mpoly_struct;

typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_zech_ctx_t fqctx;
} fq_zech_mpoly_ctx_struct;

typedef fq_zech_mpoly_ctx_struct fq_zech_mpoly_ctx_t[1];

typedef struct
{
    fq_zech_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;     /* number of bits per exponent */
} fq_zech_mpoly_struct;

typedef fq_zech_mpoly_struct fq_zech_mpoly_t[1];

typedef struct
{
    mp_limb_t constant;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} nmod_mpoly_factor_struct;

typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

typedef struct
{
    fmpz_t constant;
    fmpz_t constant_den;        /* should be one after normal operations */
    fmpz_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fmpz_mpoly_factor_struct;

typedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

typedef struct
{
    fq_nmod_t constant;
    fq_nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fq_nmod_mpoly_factor_struct;

typedef fq_nmod_mpoly_factor_struct fq_nmod_mpoly_factor_t[1];

typedef struct
{
    fmpz_mod_poly_t *polys;
    ulong p;
    ulong q;
    fmpz_mod_ctx_t ctx;
} _unity_zpq;

typedef _unity_zpq unity_zpq[1];

typedef struct
{
    fmpz_mod_poly_t poly;
    ulong p;
    ulong exp;
    fmpz_mod_ctx_t ctx;
} _unity_zp;

typedef _unity_zp unity_zp[1];

typedef enum
{
    GRAM,
    Z_BASIS
} rep_type;

typedef enum
{
    APPROX,
    EXACT
} gram_type;

typedef struct
{
    double delta;
    double eta;
    rep_type rt;
    gram_type gt;
} fmpz_lll_struct;

typedef fmpz_lll_struct fmpz_lll_t[1];

typedef struct
{
    mp_limb_t * coeffs;
    slong alloc;
    slong length;
} n_poly_struct;

typedef n_poly_struct n_poly_t[1];



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

/* The largest bit count for an fmpz to be small */
#define SMALL_FMPZ_BITCOUNT_MAX (FLINT_BITS - 2)

/* maximum positive value a small coefficient can have */
#define COEFF_MAX ((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1))

/* minimum negative value a small coefficient can have */
#define COEFF_MIN (-((WORD(1) << SMALL_FMPZ_BITCOUNT_MAX) - WORD(1)))


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
#define FLINT_ABS(x) ((slong)(x) < 0 ? -(x) : (x))
#define FLINT_SIGN_EXT(x) (-(ulong)((slong)(x) < 0))
#define FLINT_SGN(x) ((0 < (slong)(x)) - ((slong)(x) < 0))

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

/* common usage of flint_malloc */
#define FLINT_ARRAY_ALLOC(n, T) (T *) flint_malloc((n)*sizeof(T))
#define FLINT_ARRAY_REALLOC(p, n, T) (T *) flint_realloc(p, (n)*sizeof(T))

/* temporary allocation */
#define TMP_INIT                    \
    typedef struct __tmp_struct     \
    {                               \
        void * block;               \
        struct __tmp_struct * next; \
    } __tmp_t;                      \
    __tmp_t * __tmp_root;           \
    __tmp_t * __tpx

#define TMP_START                   \
    __tmp_root = NULL

#if FLINT_WANT_ASSERT
#define TMP_ALLOC(size)                             \
    (__tpx = (__tmp_t *) alloca(sizeof(__tmp_t)),   \
     __tpx->next = __tmp_root,                      \
     __tmp_root = __tpx,                            \
     __tpx->block = flint_malloc(size))
#else
#define TMP_ALLOC(size)                             \
   (((size) > 8192) ?                               \
      (__tpx = (__tmp_t *) alloca(sizeof(__tmp_t)), \
       __tpx->next = __tmp_root,                    \
       __tmp_root = __tpx,                          \
       __tpx->block = flint_malloc(size)) :         \
      alloca(size))
#endif

#define TMP_ARRAY_ALLOC(n, T) (T *) TMP_ALLOC((n)*sizeof(T))

#define TMP_END                         \
    while (__tmp_root)                  \
    {                                   \
        flint_free(__tmp_root->block);  \
        __tmp_root = __tmp_root->next;  \
    }

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
