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

/* printing */
#if defined(_WIN64) || defined(__mips64)
#if defined(__MINGW64__)
#define WORD_FMT "%I64"
#define WORD_WIDTH_FMT "%*I64"
#else
#define WORD_FMT "%ll"
#define WORD_WIDTH_FMT "%*ll"
#endif
#else
#define WORD_FMT "%l"
#define WORD_WIDTH_FMT "%*l"
#endif

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
 * misc
 */
#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)



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


/* common usage of flint_malloc */
#define FLINT_ARRAY_ALLOC(n, T) (T *) flint_malloc((n)*sizeof(T))
#define FLINT_ARRAY_REALLOC(p, n, T) (T *) flint_realloc(p, (n)*sizeof(T))


/* temporary allocation
 * NOTE: Requires user to include alloca */
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


#define FLINT_MPZ_PTR_SWAP(a, b)    \
  do {                              \
    mpz_ptr __tmp = (a);            \
    (a) = (b);                      \
    (b) = __tmp;                    \
  } while (0)

#define FLINT_MPN_ZERO(xxx, nnn)                \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < (nnn); ixxx++)    \
            (xxx)[ixxx] = UWORD(0);             \
    } while (0)

#define FLINT_MPN_COPYI(xxx, yyy, nnn)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < (nnn); ixxx++)    \
            (xxx)[ixxx] = (yyy)[ixxx];          \
    } while (0)

#define FLINT_MPN_COPYD(xxx, yyy, nnn)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = nnn - 1; ixxx >= 0; ixxx--) \
            (xxx)[ixxx] = (yyy)[ixxx];          \
    } while (0)

#define FLINT_MPN_STORE(xxx, nnn, yyy)          \
    do                                          \
    {                                           \
        slong ixxx;                             \
        for (ixxx = 0; ixxx < nnn; ixxx++)      \
        (xxx)[ixxx] = yyy;                      \
    } while (0)

#endif
