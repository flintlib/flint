/*
    Copyright (C) 2015 Nitin Kumar
    Copyright (C) 2016 William Hart
    Copyright (C) 2020 Dan Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "qsieve.h"

/* Input for which qsieve_factor would take forever -- see #2250 */
static const char *qsieve_hang_input[] = {
    "2" "100838095593318077623227274510379021",
    /* the next two examples have moved up to the front of the list
        because they hit the rare "O(s) pre-check" branch */
    "2" "1000000000000000000000000000000000000000420217",
    "4" "8312897460409825844798204414148551477617523587847",
    "2" "311998636618373748658605151670132867",
    "2" "126582278481012658227848101265822806309",
    "2" "322580645161290322580645161290322654733",
    "2" "555555555555555555555555555555555693427",
    "2" "561482313307130825379000561482313307499",
    "2" "720461095100864553314121037463976947183",
    "2" "1628664495114006514657980456026058632341",
    "2" "1650516364044491319109183307997907145251",
    "3" "1666666666666666666666666666666666845161",
    "2" "1666666666666666666666666666666666979801",
    "2" "1672240802675585284280936454849498335537",
    "2" "64874887389308019117684712840688961808594013",
    "2" "95662720772928579912788805098723834088195937",
    "2" "97349439977949650929026562827088794528412033",
    "2" "111111111111111111111111111111111111111135749",
    "2" "277315585135884636716583471991125901275651731",
    "3" "333333333333333333333333333333333333333382979",
    "2" "500000000000000000000000000000000000000017711",
    "2" "500000000000000000000000000000000000000031439",
    "2" "500000000000000000000000000000000000000085101",
    "2" "500000000000000000000000000000000000000126443",
    "2" "500000000000000000000000000000000000000178139",
    "3" "571755288736420811892510005717552887364208149",
    "3" "633713561470215462610899873257287705956907549",
    "2" "970931288163667945383173178217350736305742379",
    "3" "1000000000000000000000000000000000000000035422",
    "2" "1000000000000000000000000000000000000000042223",
    "3" "1000000000000000000000000000000000000000062878",
    "2" "1000000000000000000000000000000000000000078147",
    "4" "1000000000000000000000000000000000000000148937",
    "2" "1000000000000000000000000000000000000000169537",
    "3" "1000000000000000000000000000000000000000170202",
    "2" "1000000000000000000000000000000000000000186843",
    "3" "1000000000000000000000000000000000000000221741",
    "3" "1000000000000000000000000000000000000000252886",
    "2" "1000000000000000000000000000000000000000253483",
    "3" "1000000000000000000000000000000000000000280147",
    "3" "1000000000000000000000000000000000000000287517",
    "2" "1000000000000000000000000000000000000000322347",
    "3" "1000000000000000000000000000000000000000356278",
    "6" "1000000000000000000000000000000000000000368621",
    "2" "1000000000000000000000000000000000000000414187",
    "2" "1000000000000000000000000000000000000000415947",
    "2" "1000000000000000000000000000000000000000449343",
    "2" "1000000000000000000000000000000000000000484747",
    "2" "1000000000000000000000000000000000000000489307",
    "2" "1000000000000000000000000000000000000000496267",
    "2" "1030927835051546391752577319587628865979381499",
    "2" "1441266465979343215563261116966752565379363587",
    "2" "1537794372287715175877542358545984665114519547",
    "5" "73644707908281013635796769783360384546398094566",
    "3" "12195121951219512195121951219512195121951219512199",
    "3" "23809523809523809523809523809523809523809523811539",
    NULL,
};

void randprime(fmpz_t p, flint_rand_t state, slong bits)
{
    fmpz_randbits(p, state, bits);

    if (fmpz_sgn(p) < 0)
       fmpz_neg(p, p);

    if (fmpz_is_even(p))
       fmpz_add_ui(p, p, 1);

    while (!fmpz_is_probabprime(p))
       fmpz_add_ui(p, p, 2);
}

TEST_FUNCTION_START(qsieve_factor, state)
{
   slong i;
   fmpz_t n, x, y, z;
   fmpz_factor_t factors;
   slong max_threads = 5;
   slong tmul = 3;
#ifdef _WIN32
   tmul = 1;
#endif

   fmpz_init(x);
   fmpz_init(y);
   fmpz_init(z);
   fmpz_init(n);

   /* Test n with large prime factor */
   {
      fmpz_set_str(n, "12387192837918273918723981291837121933111751252512531193171", 10);

      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 5)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test n with large prime factor\n");
         flint_printf("%wd factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

    for (i = 0; i < 1 + flint_test_multiplier() && qsieve_hang_input[i] != NULL; i++)
    {
        fmpz_factor_init(factors);
        fmpz_set_str(n, qsieve_hang_input[i] + 1, 0);
        slong mult = qsieve_hang_input[i][0] - '0';
        /* qsieve_factor is tested indirectly */
        fmpz_factor(factors, n);
        if (mult != factors->num)
        {
             flint_printf("FAIL:\n");
             flint_printf("Table entry %{fmpz}, \n%wd factors,\n", n, mult);
             flint_printf("%wd factors found\n", factors->num);
             fflush(stdout);
             flint_abort();
        }
        fmpz_factor_clear(factors);
    }

   /* Test random n, two factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      slong bits = 40;

      randprime(x, state, bits);
      do {
         randprime(y, state, bits);
      } while (fmpz_equal(x, y));

      fmpz_mul(n, x, y);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 2)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, two factors\ni = %wd\n", i);
         flint_printf("%wd factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   /* Test random n, three factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      randprime(x, state, 40);
      do {
         randprime(y, state, 40);
      } while (fmpz_equal(x, y));
      do {
         randprime(z, state, 40);
      } while (fmpz_equal(x, z) || fmpz_equal(y, z));

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, three factors\ni = %wd\n", i);
         flint_printf("%wd factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   /* Test random n, small factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      randprime(x, state, 10);
      do {
         randprime(y, state, 10);
      } while (fmpz_equal(x, y));
      randprime(z, state, 40);

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, small factors\ni = %wd\n", i);
         flint_printf("%wd factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);
   fmpz_clear(z);

   TEST_FUNCTION_END(state);
}
