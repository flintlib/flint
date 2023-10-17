/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

#if FLINT64
/* n < 10^16 that pass base 2, 3, 7, 61 and 24251 sprp test */
mp_limb_t composites[] = {
    UWORD(669094855201), UWORD(1052516956501), UWORD(2007193456621),
    UWORD(2744715551581), UWORD(9542968210729), UWORD(17699592963781),
    UWORD(19671510288601), UWORD(24983920772821), UWORD(24984938689453),
    UWORD(29661584268781), UWORD(37473222618541), UWORD(46856248255981),
    UWORD(47922612926653), UWORD(48103703944453), UWORD(49110566041153),
    UWORD(49752242681221), UWORD(91206655032481), UWORD(91481980096033),
    UWORD(119034193492321), UWORD(123645258399601), UWORD(128928036060253),
    UWORD(137364148720147), UWORD(150753857310253), UWORD(153131886327421),
    UWORD(155216912613121), UWORD(185610214763821), UWORD(224334357392701),
    UWORD(227752294950181), UWORD(230058334559041), UWORD(304562854940401),
    UWORD(306001576998253), UWORD(335788261073821), UWORD(377133492079081),
    UWORD(379242177424951), UWORD(389970770948461), UWORD(397319638319521),
    UWORD(448114903362253), UWORD(523235160050221), UWORD(628999496281621),
    UWORD(699349238838253), UWORD(746667678235753), UWORD(790198268451301),
    UWORD(794036495175661), UWORD(823820871230281), UWORD(867739535711821),
    UWORD(1039918661294761), UWORD(1099127938585141), UWORD(1104388025338153),
    UWORD(1173374598605653), UWORD(1262797719066157), UWORD(1265872947674653),
    UWORD(1325898212229667), UWORD(1327034517143653), UWORD(1418575746675583),
    UWORD(1666122072463621), UWORD(1837400535259453), UWORD(1857422490084961),
    UWORD(1870756820971741), UWORD(1914550540480717), UWORD(2018963273468221),
    UWORD(2163829000939453), UWORD(2206020317369221), UWORD(2301037384029121),
    UWORD(2416062055125421), UWORD(2435076500074921), UWORD(2545656135020833),
    UWORD(2594428516569781), UWORD(2669983768115821), UWORD(2690937050990653),
    UWORD(2758640869506607), UWORD(2833525461416653), UWORD(2876662942007221),
    UWORD(2932155806957821), UWORD(2957010595723801), UWORD(3183606449929153),
    UWORD(3220133449185901), UWORD(3424103775720253), UWORD(3625360152399541),
    UWORD(3939300299037421), UWORD(3947917710714841), UWORD(3980273496750253),
    UWORD(4182256679324041), UWORD(4450605887818261), UWORD(4727893739521501),
    UWORD(4750350311306953), UWORD(4755334362931153), UWORD(5756440863559753),
    UWORD(5760976603475341), UWORD(5794399356078761), UWORD(5954850603819253),
    UWORD(6125544931991761), UWORD(6320931714094861), UWORD(6347593619672581),
    UWORD(6406268028524101), UWORD(6510632945054941), UWORD(6620082224794741),
    UWORD(6627325072566061), UWORD(6844056606431101), UWORD(6989404981060153),
    UWORD(7144293947609521), UWORD(7288348593229021), UWORD(7288539837129253),
    UWORD(7406102904971689), UWORD(7430233301822341), UWORD(7576425305871193),
    UWORD(7601696719033861), UWORD(7803926845356487), UWORD(7892007967006633),
    UWORD(7947797946559453), UWORD(8207000460596953), UWORD(8295064717807513),
    UWORD(8337196000698841), UWORD(8352714234009421), UWORD(8389755717406381),
    UWORD(8509654470665701), UWORD(8757647355282841), UWORD(8903933671696381),
    UWORD(8996133652295653), UWORD(9074421465661261), UWORD(9157536631454221),
    UWORD(9188353522314541)
};
#endif

TEST_FUNCTION_START(n_is_prime, state)
{
   int i, result;
   mp_limb_t d;
   mpz_t d_m;
   slong pow;
   ulong bits;

   /* Test that primes pass the test */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (mpz_size(d_m) > 1);

      result = n_is_prime(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared composite\n", d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   /* Test that composites do not pass */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         if (d == UWORD(1)) d++;
         flint_mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      result = !n_is_prime(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared prime\n", d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   /* Test that powers do not pass */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      pow = n_randint(state, 6) + 2;
      bits = n_randint(state, FLINT_BITS) + 1;
      bits /= pow;

      d = n_randbits(state, bits);
      d = n_pow(d, pow);

      result = !n_is_prime(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("Perfect power d = %wu is declared prime\n", d);
         fflush(stdout);
         flint_abort();
      }
   }

#if FLINT64
   for (i = 0; i < sizeof(composites) / sizeof(mp_limb_t); i++)
   {
      d = composites[i];

      result = !n_is_prime(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("Known composite d = %wu is declared prime\n", d);
         fflush(stdout);
         flint_abort();
      }
   }
#endif

   TEST_FUNCTION_END(state);
}
