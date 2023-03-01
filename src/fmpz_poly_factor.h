/*
    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2016 William Hart
    Copyright (C) 2009, 2011 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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

FLINT_DLL void fmpz_poly_factor_zassenhaus_recombination(
            fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac,
                               const fmpz_poly_t F, const fmpz_t P, slong exp);
    
FLINT_DLL void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, 
                                                          const fmpz_poly_t F);

FLINT_DLL void fmpz_poly_factor_mignotte(fmpz_t B, const fmpz_poly_t f);

FLINT_DLL void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, 
              slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij);

FLINT_DLL void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, 
                                                          const fmpz_poly_t G);

FLINT_DLL void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac,
                                               const fmpz_poly_t f, slong exp);

FLINT_DLL void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac,
                                               const fmpz_poly_t f, slong exp);

FLINT_DLL slong _fmpz_poly_factor_CLD_mat(fmpz_mat_t res, const fmpz_poly_t f,
                             fmpz_poly_factor_t lifted_fac, fmpz_t P, ulong k);

FLINT_DLL int fmpz_poly_factor_van_hoeij_check_if_solved(fmpz_mat_t M,
          fmpz_poly_factor_t final_fac, fmpz_poly_factor_t lifted_fac,
                          const fmpz_poly_t f, fmpz_t P, slong exp, fmpz_t lc);

FLINT_DLL void fmpz_poly_factor_van_hoeij(fmpz_poly_factor_t final_fac, 
        const nmod_poly_factor_t fac, const fmpz_poly_t f, slong exp, ulong p);

FLINT_DLL void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G);

/* Inlines *******************************************************************/

FLINT_DLL void fmpz_poly_factor_get_fmpz_poly(fmpz_poly_t z, const fmpz_poly_factor_t F, slong i);
FLINT_DLL void fmpz_poly_factor_get_fmpz(fmpz_t z, const fmpz_poly_factor_t F);

/* zassenhaus ****************************************************************/

FLINT_DLL void zassenhaus_subset_first(slong * s, slong r, slong m);

FLINT_DLL int zassenhaus_subset_next(slong * s, slong r);

FLINT_DLL slong zassenhaus_subset_next_disjoint(slong * s, slong r);

typedef struct {
    slong deg;
    unsigned char * pos_degs;   /* possible degrees: entries are 0 or 1*/
    slong new_length;
    slong new_total;
    slong * new_degs;
    slong alloc;
} zassenhaus_prune_struct;

typedef zassenhaus_prune_struct zassenhaus_prune_t[1];

FMPZ_POLY_FACTOR_INLINE
void zassenhaus_prune_init(zassenhaus_prune_t Z)
{
    Z->deg = 0;
    Z->pos_degs = NULL;
    Z->new_length = 0;
    Z->new_total = 0;
    Z->new_degs = NULL;
    Z->alloc = 0;
}

FLINT_DLL void zassenhaus_prune_clear(zassenhaus_prune_t Z);

FLINT_DLL void zassenhaus_prune_set_degree(zassenhaus_prune_t Z, slong d);

FMPZ_POLY_FACTOR_INLINE
void zassenhaus_prune_start_add_factors(zassenhaus_prune_t Z)
{
    Z->new_length = 0;
    Z->new_total = 0;
}

FLINT_DLL void zassenhaus_prune_add_factor(zassenhaus_prune_t Z,
                                                         slong deg, slong exp);

FLINT_DLL void zassenhaus_prune_end_add_factors(zassenhaus_prune_t Z);

FLINT_DLL int zassenhaus_prune_must_be_irreducible(const zassenhaus_prune_t Z);

FMPZ_POLY_FACTOR_INLINE
int zassenhaus_prune_degree_is_possible(const zassenhaus_prune_t Z, slong d)
{
    if (d <= 0)
        return d == 0;

    if (d >= Z->deg)
        return d == Z->deg;

    return Z->pos_degs[d];
}

FLINT_DLL void fmpz_poly_factor_zassenhaus_recombination_with_prune(
            fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac,
                                const fmpz_poly_t F, const fmpz_t P, slong exp,
                                                   const zassenhaus_prune_t Z);

#ifdef __cplusplus
}
#endif

#endif

