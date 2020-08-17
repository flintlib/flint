/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, 
                                         const fmpz_t m, const fmpz_t n)
{
    fmpz_t t;
    slong i, B = fmpz_sizeinbase(m, 2);

    fmpz_init(t);
    fmpz_set_ui(Vm, 2);
    fmpz_set(Vm1, A);

    for (i = B - 1; i >= 0; i--)
    {
       if (fmpz_tstbit(m, i)) /* 1 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_sub(t, t, A);
          fmpz_mod(Vm, t, n);

          fmpz_mul(t, Vm1, Vm1);
          fmpz_sub_ui(t, t, 2);
          fmpz_mod(Vm1, t, n);
       } else /* 0 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_sub(t, t, A);
          fmpz_mod(Vm1, t, n);

          fmpz_mul(t, Vm, Vm);
          fmpz_sub_ui(t, t, 2);
          fmpz_mod(Vm, t, n);
       }
    }

    fmpz_clear(t);
}

void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, 
                                         const fmpz_t m, const fmpz_t n)
{
    fmpz_t t, Q;
    slong i, bits = fmpz_sizeinbase(m, 2);

    fmpz_init(t);
    fmpz_init(Q);
    fmpz_set_ui(Q, 1);
    fmpz_set_ui(Vm, 2);
    fmpz_set(Vm1, A);

    for (i = bits - 1; i >= 0; i--)
    {
       if (fmpz_tstbit(m, i)) /* 1 in binary repn */
       {
          fmpz_mul(t, Vm1, Vm);
          fmpz_submul(t, Q, A);
          fmpz_mod(Vm, t, n);

          fmpz_mul(Vm1, Vm1, Vm1);
          fmpz_mul_ui(t, Q, 2);
          fmpz_mul(t, t, B);
          fmpz_sub(Vm1, Vm1, t);
          fmpz_mod(Vm1, Vm1, n);

          fmpz_mul(Q, Q, Q);
          fmpz_mul(Q, Q, B);
          fmpz_mod(Q, Q, n);
       } else /* 0 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_submul(t, Q, A);
          fmpz_mod(Vm1, t, n);

          fmpz_mul(t, Vm, Vm);
          fmpz_submul_ui(t, Q, 2);
          fmpz_mod(Vm, t, n);

          fmpz_mul(Q, Q, Q);
          fmpz_mod(Q, Q, n);
       }
    }

    fmpz_clear(Q);
    fmpz_clear(t);
}

/* Compute U_{2m}, U_{2m + 1} given U_m, U_{m + 1} */
void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, 
                            const fmpz_t Um1, const fmpz_t A, const fmpz_t B,
                            const fmpz_t n)
{
   fmpz_t t, t2;

   fmpz_init(t);
   fmpz_init(t2);

   fmpz_mul_2exp(t, Um1, 1); /* U_m(2U_{m+1) - AU_m) */
   fmpz_submul(t, Um, A);
   fmpz_mul(t, t, Um);

   fmpz_mul(U2m1, Um1, Um1); /* U_{m+1}^2 - BU_m^2 */
   fmpz_mul(t2, Um, Um);
   fmpz_submul(U2m1, t2, B);
   fmpz_mod(U2m1, U2m1, n);

   fmpz_mod(U2m, t, n);

   fmpz_clear(t);
   fmpz_clear(t2);
}

/* 
   Compute U_{m + n}, U_{m + n + 1} given U_m, U_{m + 1} and
   U_n, U_{n + 1}
*/
void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, 
                            const fmpz_t Um1, const fmpz_t Un, 
                            const fmpz_t Un1, const fmpz_t A, const fmpz_t B, 
                            const fmpz_t n)
{
   fmpz_t t, t2;

   fmpz_init(t);
   fmpz_init(t2);

   fmpz_mul(t, Un, A); /* U_nU_{m + 1} - BU_m(AU_n - U_{n + 1})/B */ 
   fmpz_sub(t, Un1, t);
   fmpz_mul(t, t, Um);
   fmpz_addmul(t, Un, Um1);

   fmpz_mul(Umn1, Un1, Um1); /* U_{n + 1}U_{m + 1} - BU_mU_n */
   fmpz_mul(t2, Um, Un);
   fmpz_submul(Umn1, t2, B);
   fmpz_mod(Umn1, Umn1, n);

   fmpz_mod(Umn, t, n);

   fmpz_clear(t);
   fmpz_clear(t2);
}

/* Compute U_{km}, U_{km + 1} from U_m, U_{m + 1}, k > 0 */
void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1,
                        const fmpz_t Um, const fmpz_t Um1,
                         const fmpz_t A, const fmpz_t B, const fmpz_t k, 
                         const fmpz_t n)
{
   slong i = 0, b = fmpz_sizeinbase(k, 2);
   fmpz_t t, t1;

   fmpz_init(t);
   fmpz_init(t1);

   fmpz_set(Ukm, Um);
   fmpz_set(Ukm1, Um1);

   while (!fmpz_tstbit(k, i))
   {
      fmpz_lucas_chain_double(Ukm, Ukm1, Ukm, Ukm1, A, B, n);
      i++;
   }

   i++;
   
   if (i < b)
   {
      fmpz_set(t, Ukm);
      fmpz_set(t1, Ukm1);
   }

   while (i < b)
   {
      fmpz_lucas_chain_double(t, t1, t, t1, A, B, n);
      if (fmpz_tstbit(k, i))
         fmpz_lucas_chain_add(Ukm, Ukm1, Ukm, Ukm1, t, t1, A, B, n);
      i++;
   }

   fmpz_clear(t);
   fmpz_clear(t1);
}

/* Compute U_m, U_{m + 1} from V_m, V_{m + 1} */
void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, 
                           const fmpz_t Vm, const fmpz_t Vm1,
                           const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, 
                           const fmpz_t n)
{
   fmpz_t t;

   fmpz_init(t);

   fmpz_mul_2exp(t, Vm1, 1); /* (2V_{m + 1} - AV_m) / D */
   fmpz_submul(t, Vm, A);
   fmpz_mul(t, t, Dinv);
   
   fmpz_set(Um1, Vm);
   fmpz_mod(Um, t, n);

   fmpz_addmul(Um1, Um, A); /* (V_m + AU_m) / 2 */
   if (!fmpz_is_even(Um1))
      fmpz_add(Um1, Um1, n);
   fmpz_tdiv_q_2exp(Um1, Um1, 1);

   fmpz_mod(Um1, Um1, n);

   fmpz_clear(t);
}
