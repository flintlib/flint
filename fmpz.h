/*============================================================================

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

===============================================================================*/
/******************************************************************************

 Copyright (C) 2009 William Hart
 
******************************************************************************/

#ifndef FMPZ_H
#define FMPZ_H

#include <stdio.h>
#include <mpir.h>

typedef long fmpz;
typedef fmpz fmpz_t[1];

// maximum positive value a small coefficient can have
#define COEFF_MAX ((1L<<(FLINT_BITS-2))-1L)

// minimum negative value a small coefficient can have
#define COEFF_MIN (-((1L<<(FLINT_BITS-2))-1L))

// turn a pointer to an __mpz_struct into a fmpz_t
#define PTR_TO_COEFF(xxx) ((((ulong)xxx)>>2) | (1L<<(FLINT_BITS - 2))) 

// turns an fmpz into a pointer to an mpz
#define COEFF_TO_PTR(xxx) ((__mpz_struct *) (xxx<<2))

#define COEFF_IS_MPZ(xxx) ((xxx>>(FLINT_BITS-2)) == 1L) // is xxx a pointer not an integer

extern gmp_randstate_t fmpz_randstate;

__mpz_struct * _fmpz_new_mpz(void);

void _fmpz_clear_mpz(fmpz f);

void _fmpz_cleanup(void);

__mpz_struct * _fmpz_promote(fmpz_t f);

__mpz_struct * _fmpz_promote_val(fmpz_t f);

static inline
void _fmpz_demote(fmpz_t f)
{
	if (COEFF_IS_MPZ(*f)) 
	{
		_fmpz_clear_mpz(*f);
		(*f) = 0L;
	}
}

void _fmpz_demote_val(fmpz_t f);

static inline
void fmpz_init(fmpz_t f)
{
	(*f) = 0L;
}

void fmpz_init2(fmpz_t f, ulong limbs);

static inline
void fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

void fmpz_randinit(void);

void fmpz_randclear(void);

void fmpz_randtest(fmpz_t f, mp_bitcnt_t bits);

long fmpz_get_si(const fmpz_t f);

ulong fmpz_get_ui(const fmpz_t f);

void fmpz_set_si(fmpz_t f, const long val);

void fmpz_set_ui(fmpz_t f, const ulong val);

void fmpz_get_mpz(mpz_t x, const fmpz_t f);

void fmpz_set_mpz(fmpz_t f, const mpz_t x);

static inline
void fmpz_zero(fmpz_t f)
{
   _fmpz_demote(f);	
   (*f) = 0L;
}

void fmpz_set(fmpz_t f, const fmpz_t g);

int fmpz_equal(const fmpz_t f, const fmpz_t g);

static inline 
void fmpz_print(fmpz_t x)
{
	if (!COEFF_IS_MPZ(*x)) printf("%ld", *x);
	else gmp_printf("%Zd", COEFF_TO_PTR(*x));
}

static inline
void fmpz_swap(fmpz_t f, fmpz_t g)
{
	if (f == g) return; // swapping not required
	
	fmpz t = *f;
   *f = *g;
	*g = t;
}

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);

int fmpz_cmp(const fmpz_t f, const fmpz_t g);

mp_size_t fmpz_size(const fmpz_t f);

int fmpz_sgn(const fmpz_t f);

mp_bitcnt_t fmpz_bits(const fmpz_t f);

void fmpz_neg(fmpz_t f1, const fmpz_t f2);

void fmpz_abs(fmpz_t f1, const fmpz_t f2);

void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_sub(fmpz_t f, const fmpz_t g, fmpz_t h);

void fmpz_mul_ui(fmpz_t f, const fmpz_t g, const ulong x);

void fmpz_mul_si(fmpz_t f, const fmpz_t g, const long x);

void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, const ulong exp);

void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, const ulong exp);

void fmpz_add_ui(fmpz_t f, const fmpz_t g, const ulong x);

void fmpz_sub_ui(fmpz_t f, const fmpz_t g, const ulong x);

void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, const ulong x);

void fmpz_submul_ui(fmpz_t f, const fmpz_t g, const ulong x);

#endif







