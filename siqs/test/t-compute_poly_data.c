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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

int main(void)
{
   int i, j, k;
   mp_limb_t small_factor, p, pinv, pmod, mod_inv, q, qsqrt;
   fmpz_t n;
   fmpz_t temp, temp2, B, B_2;
   qs_t qs_inf;
   fmpz_init(n);
   fmpz_init(temp);
   fmpz_init(temp2);
   fmpz_init(B_2);
   FLINT_TEST_INIT(state);

   flint_printf("compute_poly_data....");
   fflush(stdout);

   for (i = 0; i < 1000; i++)
   {
       fmpz_randtest_unsigned(n, state, 270);

       if (fmpz_is_zero(n) || fmpz_is_one(n) || fmpz_bits(n) <= 200) continue;

       qsieve_init(qs_inf, n);
       small_factor = qsieve_knuth_schroeppel(qs_inf);

       if (small_factor) continue;

       fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k); /* haven't calculated earlier */
       small_factor = qsieve_primes_init(qs_inf);

       if (small_factor) continue;

       qsieve_poly_init(qs_inf);

       small_factor = qsieve_compute_q0(qs_inf);

       if (small_factor) continue;

       qsieve_init_A0(qs_inf);

       qsieve_compute_pre_data(qs_inf);

       /* check correctness of precomputed data for 'A0' */

       /* check values of 'A0_divp' */
       for (j = 0; j < qs_inf->s; j++)
       {
           p = qs_inf->factor_base[qs_inf->A_ind[j]].p;
           fmpz_divexact_ui(temp, qs_inf->A0, p);

           if (fmpz_cmp(qs_inf->A0_divp[j], temp) != 0)
           {
               fmpz_print(qs_inf->A0);
               flint_printf(" / %wu != ", p);
               fmpz_print(qs_inf->A0_divp[j]);
               flint_printf("\n");
               abort();
           }
       }

       /* check values of 'A0_inv' */
       for (j = 2; j < qs_inf->num_primes; j++)
       {
           p = qs_inf->factor_base[j].p;
           fmpz_mul_ui(temp, qs_inf->A0, qs_inf->A0_inv[j]);
           pmod = fmpz_fdiv_ui(temp, p);

           if (pmod != UWORD(0) && pmod != UWORD(1))
           {
               fmpz_print(qs_inf->A0);
               flint_printf(" ^ (%wd) modulo %wu != %wu\n", -1, p, qs_inf->A_inv[j]);
               abort();
           }
       }

      /* check values of 'B0_terms' */
       for (j = 0; j < qs_inf->s; j++)
       {
           p = qs_inf->factor_base[qs_inf->A_ind[j]].p;
           fmpz_mul_ui(temp, qs_inf->A0_divp[j], qs_inf->B0_terms[j]);
           fmpz_pow_ui(B_2, temp, UWORD(2));

           if (fmpz_fdiv_ui(qs_inf->kn, p) != fmpz_fdiv_ui(B_2, p))
           {
               flint_printf("( ");
               fmpz_print(qs_inf->A0_divp[j]);
               flint_printf(" * %wu )", qs_inf->B0_terms[j]);
               flint_printf(" ^ %wd modulo %wu != ", UWORD(2), p);
               fmpz_print(qs_inf->kn);
               flint_printf(" modulo %wu \n", p);
               abort();
           }

           for (k = 0; k < qs_inf->s; k++)
           {
               if (k == j) continue;

               q = qs_inf->factor_base[qs_inf->A_ind[k]].p;

               if (fmpz_fdiv_ui(temp, q) != UWORD(0))
               {
                   flint_printf("( ");
                   fmpz_print(qs_inf->A0_divp[j]);
                   flint_printf(" * %wu )", qs_inf->B0_terms[j]);
                   flint_printf(" ^ %wd modulo %wu != %wu\n", UWORD(2), q, UWORD(0));
                   abort();
               }
           }
       }

       qsieve_init_poly_first(qs_inf);

       /*
          check correctness of precomputed data corresponding to initialization
          of first polynomial
       */

       /* check values of 'A_divp' */
       for (j = 0; j < qs_inf->s; j++)
       {
           p = qs_inf->factor_base[qs_inf->A_ind[j]].p;
           fmpz_divexact_ui(temp, qs_inf->A, p);

           if (fmpz_cmp(qs_inf->A_divp[j], temp) != 0)
           {
               fmpz_print(qs_inf->A);
               flint_printf(" / %wu != ", p);
               fmpz_print(qs_inf->A_divp[j]);
               flint_printf("\n");
               abort();
           }
       }

       /* check values of 'A_inv' */
       for (j = 2; j < qs_inf->num_primes; j++)
       {
           p = qs_inf->factor_base[j].p;
           fmpz_mul_ui(temp, qs_inf->A, qs_inf->A_inv[j]);
           pmod = fmpz_fdiv_ui(temp, p);

           if (pmod != UWORD(0) && pmod != UWORD(1))
           {
               fmpz_print(qs_inf->A);
               flint_printf(" ^ (%wd) modulo %wu != %wu\n", -1, p, qs_inf->A_inv[j]);
               abort();
           }
       }

       /* check values of 'B_terms' */
       for (j = 0; j < qs_inf->s; j++)
       {
           p = qs_inf->factor_base[qs_inf->A_ind[j]].p;
           fmpz_pow_ui(temp, qs_inf->B_terms[j], UWORD(2));

           if (fmpz_fdiv_ui(qs_inf->kn, p) != fmpz_fdiv_ui(temp, p))
           {
               fmpz_print(temp);
               flint_printf(" ^ %wd modulo %wu != ", UWORD(2), p);
               fmpz_print(qs_inf->kn);
               flint_printf(" modulo %wu \n", p);
               abort();
           }

           for (k = 0; k < qs_inf->s; k++)
           {
               if (k == j) continue;

               q = qs_inf->factor_base[qs_inf->A_ind[k]].p;

               if (fmpz_fdiv_ui(temp, q) != UWORD(0))
               {
                   flint_printf("( ");
                   fmpz_print(qs_inf->A_divp[j]);
                   flint_printf(" * %wu )", qs_inf->B_terms[j]);
                   flint_printf(" ^ %wd modulo %wu != %wu\n", UWORD(2), q, UWORD(0));
                   abort();
               }
           }
       }

       /* check for initial 'B' coefficient */

       fmpz_pow_ui(B_2, qs_inf->B[0], UWORD(2));
       fmpz_mod(temp, B_2, qs_inf->A);
       fmpz_mod(temp2, qs_inf->kn, qs_inf->A);

       if (fmpz_cmp(temp, temp2) != 0)
       {
           fmpz_print(qs_inf->B[0]);
           flint_printf(" ^ %wu modulo ", UWORD(2));
           fmpz_print(qs_inf->A);
           flint_printf(" != ");
           fmpz_print(qs_inf->kn);
           flint_printf(" modulo ");
           fmpz_print(qs_inf->A);
           flint_printf("\n");
           abort();
       }

       /* check values of 'A_inv2B' */
       for (j =0; j < qs_inf->s; j++)
       {
           for (k = 2; k < qs_inf->num_primes; k++)
           {
               p = qs_inf->factor_base[k].p;
               pinv = qs_inf->factor_base[k].pinv;
               pmod = fmpz_fdiv_ui(qs_inf->B_terms[j], p);
               pmod *= 2;
               pmod = n_mulmod2_preinv(pmod, qs_inf->A_inv[k], p, pinv);

               if (pmod != qs_inf->A_inv2B[j][k])
               {
                   flint_printf("( %wu * ",UWORD(2));
                   fmpz_print(qs_inf->B_terms[j]);
                   flint_printf(" * %wu ) modulo %wu != %wu",
                                qs_inf->A_inv[k], p, qs_inf->A_inv2B[j][k]);
                   abort();
               }
           }
       }

       /* check values of 'soln1' and 'soln2' */
       for (j = 2; j < qs_inf->num_primes; j++)
       {
           p = qs_inf->factor_base[j].p;
           fmpz_mul_ui(temp, qs_inf->A, qs_inf->soln1[j]);
           fmpz_add(temp2, temp, qs_inf->B[0]);
           fmpz_pow_ui(B_2, temp2, UWORD(2));

           if (fmpz_fdiv_ui(B_2, p) != fmpz_fdiv_ui(qs_inf->kn, p))
           {
               flint_printf("(");
               fmpz_print(qs_inf->A);
               flint_printf(" * %wu + ", qs_inf->soln1[j]);
               fmpz_print(qs_inf->B[0]);
               flint_printf(") ^ %wu modulo %wu != ", UWORD(2), p);
               fmpz_print(qs_inf->kn);
               flint_printf(" modulo %wu\n", p);
               abort();
           }

           fmpz_mul_ui(temp, qs_inf->A, qs_inf->soln2[j]);
           fmpz_add(temp2, temp, qs_inf->B[0]);
           fmpz_pow_ui(B_2, temp2, UWORD(2));

           if (fmpz_fdiv_ui(B_2, p) != fmpz_fdiv_ui(qs_inf->kn, p))
           {
               flint_printf("(");
               fmpz_print(qs_inf->A);
               flint_printf(" * %wu + ", qs_inf->soln2[j]);
               fmpz_print(qs_inf->B[0]);
               flint_printf(") ^ %wu modulo %wu != ", UWORD(2), p);
               fmpz_print(qs_inf->kn);
               flint_printf(" modulo %wu\n", p);
               abort();
           }
        }

        /* check for all subsequent polynomial for current 'A' */

        for (j = 1; j < (1 << qs_inf->s); j++)
        {
            qsieve_init_poly_next(qs_inf);

            /* check for current 'B' coefficient */

            fmpz_pow_ui(B_2, qs_inf->B[j], UWORD(2));
            fmpz_mod(temp, B_2, qs_inf->A);
            fmpz_mod(temp2, qs_inf->kn, qs_inf->A);

            if (fmpz_cmp(temp, temp2) != 0)
            {
                fmpz_print(qs_inf->B[j]);
                flint_printf(" ^ %wu modulo ", UWORD(2));
                fmpz_print(qs_inf->A);
                flint_printf(" != ");
                fmpz_print(qs_inf->kn);
                flint_printf(" modulo ");
                fmpz_print(qs_inf->A);
                flint_printf("\n");
                abort();
            }

            /* check for roots corresponding to current 'A' and 'B' */

            for (k = 2; k < qs_inf->num_primes; k++)
            {
                p = qs_inf->factor_base[k].p;
                fmpz_mul_ui(temp, qs_inf->A, qs_inf->soln1[k]);
                fmpz_add(temp2, temp, qs_inf->B[j]);
                fmpz_pow_ui(B_2, temp2, UWORD(2));

                if (fmpz_fdiv_ui(B_2, p) != fmpz_fdiv_ui(qs_inf->kn, p))
                {
                    flint_printf("(");
                    fmpz_print(qs_inf->A);
                    flint_printf(" * %wu + ", qs_inf->soln1[k]);
                    fmpz_print(qs_inf->B[j]);
                    flint_printf(") ^ %wu modulo %wu != ", UWORD(2), p);
                    fmpz_print(qs_inf->kn);
                    flint_printf(" modulo %wu\n", p);
                    abort();
                }

                fmpz_mul_ui(temp, qs_inf->A, qs_inf->soln2[k]);
                fmpz_add(temp2, temp, qs_inf->B[j]);
                fmpz_pow_ui(B_2, temp2, UWORD(2));

                if (fmpz_fdiv_ui(B_2, p) != fmpz_fdiv_ui(qs_inf->kn, p))
                {
                    flint_printf("(");
                    fmpz_print(qs_inf->A);
                    flint_printf(" * %wu + ", qs_inf->soln2[k]);
                    fmpz_print(qs_inf->B[j]);
                    flint_printf(") ^ %wu modulo %wu != ", UWORD(2), p);
                    fmpz_print(qs_inf->kn);
                    flint_printf(" modulo %wu\n", p);
                    abort();
                }
            }

        }

        qsieve_clear(qs_inf);
   }

   fmpz_clear(n);
   fmpz_clear(temp);
   fmpz_clear(temp2);
   fmpz_clear(B_2);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}

