/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#ifndef FMPZ_MOD_POLY_FACTOR_H
#define FMPZ_MOD_POLY_FACTOR_H

#ifdef FMPZ_MOD_POLY_FACTOR_INLINES_C
#define FMPZ_MOD_POLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_MOD_POLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Factoring  ****************************************************************/

typedef struct
{
    fmpz_mod_poly_struct *poly;
    slong *exp;
    slong num;
    slong alloc;
} fmpz_mod_poly_factor_struct;

typedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1];

FLINT_DLL void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac);

typedef struct
{
    fmpz_mod_poly_struct * baby;
    fmpz_mod_poly_struct res;
    fmpz_mod_poly_struct H;
    fmpz_mod_poly_struct v;
    fmpz_mod_poly_struct vinv;
    slong m;
}
fmpz_mod_poly_interval_poly_arg_t;

FLINT_DLL void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac);

FLINT_DLL void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac, slong alloc);

FLINT_DLL void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac, slong len);

FLINT_DLL void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac);

FLINT_DLL void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac,
                                 const fmpz_mod_poly_t poly, slong exp);

FLINT_DLL void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac);

FLINT_DLL void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res,
                                 const fmpz_mod_poly_factor_t fac);

FLINT_DLL void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, slong exp);

FLINT_DLL int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f);

FLINT_DLL int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f);

FLINT_DLL int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f);

FLINT_DLL int fmpz_mod_poly_is_irreducible_rabin_f(fmpz_t fac, 
                                                   const fmpz_mod_poly_t f);

FLINT_DLL int _fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_t p);

FLINT_DLL int _fmpz_mod_poly_is_squarefree_f(fmpz_t fac, 
                                    const fmpz * f, slong len, const fmpz_t p);

FLINT_DLL int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f);

FLINT_DLL int fmpz_mod_poly_is_squarefree_f(fmpz_t fac, const fmpz_mod_poly_t f);

FLINT_DLL int fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor,
          flint_rand_t state, const fmpz_mod_poly_t pol, slong d);

FLINT_DLL void fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors,
                               const fmpz_mod_poly_t pol, slong d);

FLINT_DLL void fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res,
                              const fmpz_mod_poly_t poly, slong * const *degs);

FLINT_DLL void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                                      const fmpz_mod_poly_t f);

FLINT_DLL void fmpz_mod_poly_factor_distinct_deg_threaded(fmpz_mod_poly_factor_t res,
                               const fmpz_mod_poly_t poly, slong * const *degs);

FLINT_DLL void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                                      const fmpz_mod_poly_t f);

FLINT_DLL void fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res,
                          const fmpz_mod_poly_t f);

FLINT_DLL void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                                      const fmpz_mod_poly_t f);

FLINT_DLL void fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                                          const fmpz_mod_poly_t poly);

FLINT_DLL void fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
                                     const fmpz_mod_poly_t f);

FLINT_DLL void * _fmpz_mod_poly_interval_poly_worker(void* arg_ptr);

#ifdef __cplusplus
}
#endif

#endif
