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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef FMPQ_SPECIAL_H
#define FMPQ_SPECIAL_H

#ifdef __cplusplus
 extern "C" {
#endif

void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k);
void fmpq_dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k);
void fmpq_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k);
void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k);

#ifdef __cplusplus
}
#endif

#endif

