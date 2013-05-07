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

#undef ulong /* interferes with system includes */
#include <stdio.h>
#define ulong unsigned long

#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Factoring  ****************************************************************/

typedef struct
{
    fmpz_mod_poly_struct *poly;
    long *exp;
    long num;
    long alloc;
} fmpz_mod_poly_factor_struct;

typedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1];

void fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac);

void fmpz_mod_poly_factor_clear(fmpz_mod_poly_factor_t fac);

void fmpz_mod_poly_factor_realloc(fmpz_mod_poly_factor_t fac, long alloc);

void fmpz_mod_poly_factor_fit_length(fmpz_mod_poly_factor_t fac, long len);

void fmpz_mod_poly_factor_set(fmpz_mod_poly_factor_t res, const fmpz_mod_poly_factor_t fac);

void fmpz_mod_poly_factor_insert(fmpz_mod_poly_factor_t fac,
                                 const fmpz_mod_poly_t poly, long exp);

void fmpz_mod_poly_factor_print(const fmpz_mod_poly_factor_t fac);

void fmpz_mod_poly_factor_concat(fmpz_mod_poly_factor_t res,
                                 const fmpz_mod_poly_factor_t fac);

void fmpz_mod_poly_factor_pow(fmpz_mod_poly_factor_t fac, long exp);

int fmpz_mod_poly_is_irreducible(const fmpz_mod_poly_t f);

int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f);

int fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f);

int
_fmpz_mod_poly_is_squarefree(const fmpz * f, long len, const fmpz_t p);

int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f);

int
fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor,
          flint_rand_t state, const fmpz_mod_poly_t pol, long d);

void
fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors,
                               const fmpz_mod_poly_t pol, long d);

void fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res,
                                       const fmpz_mod_poly_t poly, long **degs);

void fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                                      const fmpz_mod_poly_t f);

void fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res,
                          const fmpz_mod_poly_t f);

void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                                      const fmpz_mod_poly_t f);

void fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                                          const fmpz_mod_poly_t poly);

void fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
                                     const fmpz_mod_poly_t f);

#ifdef __cplusplus
}
#endif

#endif
