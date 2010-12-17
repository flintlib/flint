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

    Copyright (C) 2009 William Hart
 
******************************************************************************/

#ifndef FMPZ_H
#define FMPZ_H

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>

#include "flint.h"
#include "nmod_vec.h"

typedef long fmpz;
typedef fmpz fmpz_t[1];

typedef gmp_randstate_t fmpz_randstate_t;

extern __mpz_struct * fmpz_arr;
extern gmp_randstate_t fmpz_randstate;

/* maximum positive value a small coefficient can have */
#define COEFF_MAX ((1L << (FLINT_BITS - 2)) - 1L)

/* minimum negative value a small coefficient can have */
#define COEFF_MIN (-((1L << (FLINT_BITS - 2)) - 1L))

#if FLINT_REENTRANT

/* turn a pointer to an __mpz_struct into a fmpz_t */
#define PTR_TO_COEFF(x) (((ulong) (x) >> 2) | (1L << (FLINT_BITS - 2)))

/* turns an fmpz into a pointer to an mpz */
#define COEFF_TO_PTR(x) ((__mpz_struct *) ((x) << 2))

#else

/* turn a pointer to an __mpz_struct into a fmpz_t */
#define PTR_TO_COEFF(x) ((ulong) ((x) - fmpz_arr) | (1L << (FLINT_BITS - 2))) 

/* turns an fmpz into a pointer to an mpz */
#define COEFF_TO_PTR(x) ((__mpz_struct *) (((x) ^ (1L << (FLINT_BITS - 2))) + fmpz_arr))

#endif  /* FLINT_REENTRANT */

#define COEFF_IS_MPZ(x) (((x) >> (FLINT_BITS - 2)) == 1L)  /* is x a pointer not an integer */

__mpz_struct * _fmpz_new_mpz(void);

void _fmpz_clear_mpz(fmpz f);

void _fmpz_cleanup(void);

__mpz_struct * _fmpz_promote(fmpz_t f);

__mpz_struct * _fmpz_promote_val(fmpz_t f);

static __inline__
void _fmpz_demote(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
        (*f) = 0L;
	}
}

void _fmpz_demote_val(fmpz_t f);

static __inline__
void fmpz_init(fmpz_t f)
{
	(*f) = 0L;
}

void fmpz_init2(fmpz_t f, ulong limbs);

static __inline__
void fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

void fmpz_randinit(fmpz_randstate_t state);

void fmpz_randclear(fmpz_randstate_t state);

void fmpz_randbits(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_randm(fmpz_t f, fmpz_randstate_t state, fmpz_t m);

void fmpz_randtest(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_randtest_unsigned(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_randtest_not_zero(fmpz_t f, fmpz_randstate_t state, mp_bitcnt_t bits);

long fmpz_get_si(const fmpz_t f);

ulong fmpz_get_ui(const fmpz_t f);

void fmpz_set_si(fmpz_t f, long val);

void fmpz_set_ui(fmpz_t f, ulong val);

void fmpz_get_mpz(mpz_t x, const fmpz_t f);

void fmpz_set_mpz(fmpz_t f, const mpz_t x);

int fmpz_set_str(fmpz_t f, char * str, int b);

int fmpz_abs_fits_ui(const fmpz_t f);

static __inline__
void fmpz_zero(fmpz_t f)
{
   _fmpz_demote(f);	
   (*f) = 0L;
}

static __inline__
int fmpz_is_zero(const fmpz_t f)
{
   return (*f == 0);
}

static __inline__
int fmpz_is_one(const fmpz_t f)
{
   return (*f == 1);
}

void fmpz_set(fmpz_t f, const fmpz_t g);

int fmpz_equal(const fmpz_t f, const fmpz_t g);

void fmpz_read(fmpz_t f);

void fmpz_fread(FILE * file, fmpz_t f);

void fmpz_print(const fmpz_t x);

void fmpz_fprint(FILE * file, const fmpz_t x);

size_t fmpz_sizeinbase(const fmpz_t f, int b);

char * fmpz_get_str(char * str, int b, const fmpz_t f);

static __inline__
void fmpz_swap(fmpz_t f, fmpz_t g)
{
    if (f != g)  /* swapping required */
    {
        fmpz t = *f;
        *f = *g;
        *g = t;
    }
}

int fmpz_cmp(const fmpz_t f, const fmpz_t g);

int fmpz_cmp_ui(const fmpz_t f, ulong g);

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);

mp_size_t fmpz_size(const fmpz_t f);

int fmpz_sgn(const fmpz_t f);

mp_bitcnt_t fmpz_bits(const fmpz_t f);

void fmpz_neg(fmpz_t f1, const fmpz_t f2);

void fmpz_abs(fmpz_t f1, const fmpz_t f2);

void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_mul_si(fmpz_t f, const fmpz_t g, long x);

void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong x);

void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);

void fmpz_sqrt(fmpz_t f, const fmpz_t g);

void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

ulong fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);

int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_divexact_si(fmpz_t f, const fmpz_t g, long h);

void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);

void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);

void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h);

double fmpz_get_d_2exp(long * exp, const fmpz_t f);

int fmpz_bit_pack(mp_ptr arr, mp_bitcnt_t shift, mp_bitcnt_t bits, 
                  const fmpz_t coeff, int negate, int borrow);

int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, mp_bitcnt_t shift, 
                    mp_bitcnt_t bits, int negate, int borrow);

void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, 
                              mp_bitcnt_t shift, mp_bitcnt_t bits);

static __inline__
void fmpz_CRT_ui_precomp(fmpz_t out, fmpz_t r1, fmpz_t m1, 
                             ulong r2, ulong m2, ulong c, double pre)
{  
   ulong r1mod, s;
   fmpz_t r1modd, sm1; 
   
   fmpz_init(sm1);
   fmpz_init(r1modd);
   
   r1mod = fmpz_mod_ui(r1modd, r1, m2); 
   s = n_submod(r2, r1mod, m2);
   s = n_mulmod_precomp(s, c, m2, pre);
   
   fmpz_mul_ui(sm1, m1, s);
   fmpz_add(out, r1, sm1);
   
   fmpz_clear(sm1);
   fmpz_clear(r1modd);
}

static __inline__
void fmpz_CRT_ui2_precomp(fmpz_t out, fmpz_t r1, fmpz_t m1, 
                              ulong r2, ulong m2, ulong c, double pre)
{
   ulong r1mod, s;
   fmpz_t r1modd, sm1;

   fmpz_init(sm1);
   fmpz_init(r1modd);
   
   r1mod = fmpz_mod_ui(r1modd, r1, m2); 
   s = n_submod(r2, r1mod, m2);
   s = n_mulmod2_preinv(s, c, m2, pre); 
   
   fmpz_mul_ui(sm1, m1, s);
   fmpz_add(out, r1, sm1);

   fmpz_clear(sm1);
   fmpz_clear(r1modd);
}

void fmpz_fac_ui(fmpz_t f, ulong n);


#define FLINT_FMPZ_LOG_MULTI_MOD_CUTOFF 2

typedef struct
{
    mp_limb_t * primes;
    long num_primes;
    long n;         /* we have 2^n >= num_primes > 2^(n-1) */
    fmpz ** comb;   /* Array of arrays of products */
    fmpz ** res;    /* successive residues r_i^-1 mod r_{i+1} for pairs r_i, r_{i+1} */
    nmod_t * mod;
}
fmpz_comb_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];

fmpz ** fmpz_comb_temp_init(fmpz_comb_t comb);

void fmpz_comb_temp_free(fmpz_comb_t comb, fmpz ** comb_temp);

void fmpz_comb_init(fmpz_comb_t comb, mp_limb_t * primes, long num_primes);

void fmpz_comb_clear(fmpz_comb_t comb);

void fmpz_multi_mod_ui(mp_limb_t * out, fmpz_t in, fmpz_comb_t comb,
    fmpz ** comb_temp, fmpz_t temp);

void fmpz_multi_CRT_ui_unsigned(fmpz_t output, mp_limb_t * residues,
    fmpz_comb_t comb, fmpz ** comb_temp, fmpz_t temp, fmpz_t temp2);

void fmpz_multi_CRT_ui(fmpz_t output, mp_limb_t * residues,
    fmpz_comb_t comb, fmpz ** comb_temp, fmpz_t temp, fmpz_t temp2);

#endif

