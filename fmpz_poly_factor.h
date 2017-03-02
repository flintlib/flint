/*
    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2016 William Hart
    Copyright (C) 2009, 2011 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_FACTOR_H
#define FMPZ_POLY_FACTOR_H

#ifdef FMPZ_POLY_FACTOR_INLINES_C
#define FMPZ_POLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_POLY_FACTOR_INLINE static __inline__
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
#include "fmpz_mat.h"
#include "nmod_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

FLINT_DLL void fmpz_poly_factor_init(fmpz_poly_factor_t fac);

FLINT_DLL void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, slong alloc);

FLINT_DLL void fmpz_poly_factor_realloc(fmpz_poly_factor_t fac, slong alloc);

FLINT_DLL void fmpz_poly_factor_fit_length(fmpz_poly_factor_t fac, slong len);

FLINT_DLL void fmpz_poly_factor_clear(fmpz_poly_factor_t fac);

FLINT_DLL void fmpz_poly_factor_set(fmpz_poly_factor_t res,
                                                 const fmpz_poly_factor_t fac);

FLINT_DLL void fmpz_poly_factor_insert(fmpz_poly_factor_t fac, 
                             const fmpz_poly_t p, slong exp);

FLINT_DLL void fmpz_poly_factor_concat(fmpz_poly_factor_t res, 
                             const fmpz_poly_factor_t fac);

FLINT_DLL void fmpz_poly_factor_print(const fmpz_poly_factor_t fac);

FLINT_DLL void fmpz_poly_factor_zassenhaus_recombination(fmpz_poly_factor_t final_fac, 
	const fmpz_poly_factor_t lifted_fac, 
                               const fmpz_poly_t F, const fmpz_t P, slong exp);
    
FLINT_DLL void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, 
                                                          const fmpz_poly_t F);

FLINT_DLL void fmpz_poly_factor_mignotte(fmpz_t B, const fmpz_poly_t f);

FLINT_DLL void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, 
              slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij);

FLINT_DLL void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, 
                                                          const fmpz_poly_t G);

FLINT_DLL slong _fmpz_poly_factor_CLD_mat(fmpz_mat_t res, const fmpz_poly_t f,
                             fmpz_poly_factor_t lifted_fac, fmpz_t P, ulong k);

FLINT_DLL int fmpz_poly_factor_van_hoeij_check_if_solved(fmpz_mat_t M,
          fmpz_poly_factor_t final_fac, fmpz_poly_factor_t lifted_fac,
                          const fmpz_poly_t f, fmpz_t P, slong exp, fmpz_t lc);

FLINT_DLL void fmpz_poly_factor_van_hoeij(fmpz_poly_factor_t final_fac, 
        const nmod_poly_factor_t fac, const fmpz_poly_t f, slong exp, ulong p);

FLINT_DLL void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G);

#ifdef __cplusplus
}
#endif

#endif

