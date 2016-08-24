/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_H
#define FMPZ_MPOLY_H

#ifdef FMPZ_MPOLY_INLINES_C
#define FMPZ_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef enum {
   ORD_LEX, ORD_REVLEX, ORD_DEGLEX, ORD_DEGREVLEX
} ordering_t;

typedef struct
{
   slong N;        /* number of words per exponent vector */
   slong bits;     /* number of bits per exponent */
   slong n;        /* number of elements in exponent vector (including deg) */
   ordering_t ord; /* polynomial ordering */
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

typedef struct
{
    fmpz * coeffs; /* alloc fmpzs */
    ulong * exps;  /* N*alloc words */
    slong alloc;
    slong length;
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

/* Macros ********************************************************************/

#define degrev_from_ord(deg, rev, ord)                                \
   (deg) = (rev) = 0;                                                 \
   do {                                                               \
      switch (ord)                                                    \
      {                                                               \
      case ORD_LEX:                                                   \
         break;                                                       \
      case ORD_DEGLEX:                                                \
         (deg) = 1; break;                                            \
      case ORD_REVLEX:                                                \
         (rev) = 1; break;                                            \
      case ORD_DEGREVLEX:                                             \
         (deg) = (rev) = 1; break;                                    \
      default:                                                        \
         flint_throw(FLINT_ERROR, "Invalid ordering in fmpz_mpoly");  \
      }                                                               \
   } while (0)

/* Context object ************************************************************/

FLINT_DLL void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, 
                             slong N, slong bits, slong nvars, ordering_t ord);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)
{
   /* nothing to be done at the moment */
}

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_mpoly_init(fmpz_mpoly_t poly, fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init2(fmpz_mpoly_t poly, slong alloc, 
                                                         fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_realloc(fmpz_mpoly_t poly, slong alloc, 
                                                         fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, 
                                                         fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_clear(fmpz_mpoly_t poly, fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_normalise(fmpz_mpoly_t poly, fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_set_length(fmpz_mpoly_t poly, slong newlen, 
                                                         fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_truncate(fmpz_mpoly_t poly, slong newlen, 
                                                         fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;

        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);

        poly->length = newlen;

        _fmpz_mpoly_normalise(poly, ctx);
    }  
}

/*  Basic manipulation *******************************************************/

FLINT_DLL ulong * _fmpz_mpoly_max_degrees1(ulong * exps, slong len, 
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL ulong * fmpz_mpoly_max_degrees(fmpz_mpoly_t poly, 
                                                         fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_gen1(fmpz * poly, ulong * exps, slong i,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i,
                                                         fmpz_mpoly_ctx_t ctx);

/* Input/output **************************************************************/

FLINT_DLL void _exp_get_degrees1(ulong * expvec, ulong v, slong bits,
                                                    slong n, int deg, int rev);

FLINT_DLL char * _fmpz_mpoly_get_str_pretty1(fmpz * poly, ulong * exps, 
                  slong len, char ** x, slong bits, slong n, int deg, int rev);

FLINT_DLL char * fmpz_mpoly_get_str_pretty(fmpz_mpoly_t poly,
                                            char ** x, fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_fprint_pretty1(FILE * file, fmpz * poly, 
    ulong * exps, slong len, char ** x, slong bits, slong n, int deg, int rev);

FLINT_DLL int fmpz_mpoly_fprint_pretty(FILE * file, 
                           fmpz_mpoly_t poly, char ** x, fmpz_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
