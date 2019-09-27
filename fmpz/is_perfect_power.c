/* mpz_perfect_power_p(arg) -- Return non-zero if ARG is a perfect power,
   zero otherwise.

Copyright 1998, 1999, 2000, 2001, 2005 Free Software Foundation, Inc.

Copyright 2008 Jason Moxham
Copyright 2017 William Hart

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

/*
  We are to determine if c is a perfect power, c = a ^ b.
  Assume c is divisible by 2^n and that codd = c/2^n is odd.
  Assume a is divisible by 2^m and that aodd = a/2^m is odd.
  It is always true that m divides n.

  * If n is prime, either 1) a is 2*aodd and b = n
		       or 2) a = c and b = 1.
    So for n prime, we readily have a solution.
  * If n is factorable into the non-trivial factors p1,p2,...
    Since m divides n, m has a subset of n's factors and b = n / m.
*/

/* This is a naive approach to recognizing perfect powers.
   Many things can be improved.  In particular, we should use p-adic
   arithmetic for computing possible roots. 
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h> /* for NULL */
#undef ulong
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

static const unsigned short primes[] =
{  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
  59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,
 137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
 227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
 313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,
 419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,
 509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,
 617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,
 727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
 829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,
 947,953,967,971,977,983,991,997,0
};
#define SMALLEST_OMITTED_PRIME 1009

int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f)
{
   ulong prime, n, n2, rem;
   mpz_t u2, q;
   int exact, exp2, sgn = fmpz_sgn(f);
   slong i, uns, usize = fmpz_size(f);
   __mpz_struct * u, * r;

   if (usize == 0)
   {
      fmpz_zero(root);
      return 2;			/* consider 0 a perfect square */
   }

   if (usize == 1)
   {
      ulong r = 0;
      ulong r2;
      ulong n = fmpz_get_ui(f);
 
      int exp = n_is_perfect_power(&r, n);

      /* get higest exponent */
      while (r > 1 && (exp2 = n_is_perfect_power(&r2, r)) != 0)
      {
         exp *= exp2;
         r = r2;
      }

      if (exp == 0)
         return 0;
      else if (sgn < 0 && (exp & 1) == 0)
      {
         while ((exp & 1) == 0)
         {
            r = r*r;
            exp >>= 1;
         }
         if (exp == 1 && n != 1)
            return 0;
         else
         {
            fmpz_set_si(root, -r);
            return exp;
         }
      } else
      {
         fmpz_set_ui(root, r);
         if (sgn < 0)
            fmpz_neg(root, root);
         return exp;
      }
   }

   u = COEFF_TO_PTR(*f);
   usize = u->_mp_size;

   n2 = mpz_scan1(u, 0);

   if (n2 == 1)
      return 0;			/* 2 divides exactly once */

   if (n2 > 1 && (n2 & (n2 - 1)) == 0 && usize < 0)
      return 0;			/* 2 has power of two  multiplicity with negative u */

   uns = FLINT_ABS(usize) - n2/FLINT_BITS;
   mpz_init2(q, FLINT_BITS*uns);
   mpz_init2(u2, FLINT_BITS*uns);

   mpz_tdiv_q_2exp(u2, u, n2);

   if (n_is_prime(n2))
      goto n2prime;

   for (i = 1; primes[i] != 0; i++)
   {
      prime = primes[i];

      if (mpz_divisible_ui_p(u2, prime))	/* divisible by this prime? */
      {
         rem = mpz_tdiv_q_ui(q, u2, prime*prime);

	 if (rem != 0)
	 {
	    mpz_clear(q);
            mpz_clear(u2);
	    return 0;		/* prime divides exactly once, reject */
	 }

	 mpz_swap(q, u2);
	 for (n = 2; ; )
	 {
	    rem = mpz_tdiv_q_ui(q, u2, prime);
	    if (rem != 0)
		break;
	    mpz_swap (q, u2);
	    n++;
	 }

	 if ((n & (n - 1)) == 0 && usize < 0)
	 {
	     mpz_clear(q);
             mpz_clear(u2);
	     return 0;		/* power of two multiplicity with negative U, reject */
	 }

	 n2 = n_gcd(n2, n);
	 if (n2 == 1)
	 {
	    mpz_clear(q);
            mpz_clear(u2);
	    return 0;		/* we have multiplicity 1 of some factor */
	 }

	 if (mpz_cmpabs_ui(u2, 1) == 0)
	 {
	    mpz_clear(q);
            mpz_clear(u2);
	    
            if (usize < 0)
            {
               if ((n2 & (n2 - 1)) == 0)
                  return 0;        /* factoring completed; not consistent power */
                
               while ((n2 & 1) == 0)
                  n2 >>= 1;
            }
	    
            r = _fmpz_promote(root);
            mpz_root(r, u, n2);
            _fmpz_demote_val(root);
	    return n2;        /* factoring completed; consistent power */
	 }
         
         /* 
            as soon as n2 becomes a prime number, stop factoring
	    either we have u=x^n2 or u is not a perfect power
         */
         if (n_is_prime(n2))
	    goto n2prime;
      }
   }

   if (n2 == 0)
   {
      /* we found no factors above; have to check all values of n */
      ulong nth;

      for (nth = usize < 0 ? 3 : 2; ; nth++)
      {
         if (!n_is_prime(nth))
	    continue;

	  exact = mpz_root(q, u2, nth);

	  if (exact)
	  {
             r = _fmpz_promote(root);
             mpz_set(r, q);
             _fmpz_demote_val(root);
	     mpz_clear(q);
             mpz_clear(u2);
             return nth;
	  }

	  if (mpz_cmpabs_ui(q, SMALLEST_OMITTED_PRIME) < 0)
	  {
	     mpz_clear(q);
             mpz_clear(u2);
	     return 0;
          }
       }
    } else
    {
       unsigned long int nth;

       /* 
          we found some factors above and we just need to consider values of n
	  that divide n2
       */

       for (nth = usize < 0 ? 3 : 2; nth <= n2; nth++)
       {
          if (!n_is_prime(nth))
	     continue;

	  if (n2 % nth != 0)
	     continue;

	  exact = mpz_root(q, u, nth);

	  if (exact)
	  {
	     r = _fmpz_promote(root);
             mpz_set(r, q);
             _fmpz_demote_val(root);
	     mpz_clear(q);
             mpz_clear(u2);
	     return nth;
	  }
	  
          if (mpz_cmpabs_ui(q, SMALLEST_OMITTED_PRIME) < 0)
	  {
	     mpz_clear(q);
             mpz_clear(u2);
	     return 0;
	  }
      }

      mpz_clear(q);
      mpz_clear(u2);
      return 0;
   }
   
n2prime:

   if (n2 == 2 && usize < 0)
   {
      mpz_clear(q);
      mpz_clear(u2);
      return 0;
   }
   
   exact = mpz_root(q, u, n2);
 
   if (exact)
   {
      r = _fmpz_promote(root);
      mpz_set(r, q);
      _fmpz_demote_val(root);
      mpz_clear(q);
      mpz_clear(u2);
      return n2;
   } else
   {
      mpz_clear(q);
      mpz_clear(u2);
      return 0;
   }
}

