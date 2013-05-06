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

    Copyright (C) 2009, 2012 William Hart
   
******************************************************************************/

#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifndef FLINT64          
mp_limb_t flint_fmpz_pseudosquares[][3] = 
{
   { 17, 0, 0 },
   { 73, 0, 0 },
   { 241, 0, 0 },
   { 1009, 0, 0 },
   { 2641, 0, 0 },
   { 8089, 0, 0 },
   { 18001, 0, 0 },
   { 53881, 0, 0 },
   { 87481, 0, 0 },
   { 117049, 0, 0 },
   { 515761, 0, 0 },
   { 1083289, 0, 0 },
   { 3206641, 0, 0 },
   { 3818929, 0, 0 },
   { 9257329, 0, 0 },
   { 22000801, 0, 0 },
   { 48473881, 0, 0 },
   { 48473881, 0, 0 },
   { 175244281, 0, 0 },
   { 427733329, 0, 0 },
   { 427733329, 0, 0 },
   { 898716289u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 1720328849u, 2u, 0 },
   { 2141495009u, 5u, 0 },
   { 3553231785u, 19u, 0 },
   { 3553231785u, 19u, 0 },
   { 2991566689u, 45u, 0 },
   { 2991566689u, 45u, 0 },
   { 2804689073u, 668u, 0 },
   { 2804689073u, 668u, 0 },
   { 2804689073u, 668u, 0 },
   { 46910577u, 6112u, 0 },
   { 46910577u, 6112u, 0 },
   { 1079027281u, 26178u, 0 },
   { 1079027281u, 26178u, 0 },
   { 1079027281u, 26178u, 0 },
   { 3590018425u, 41661u, 0 },
   { 3590018425u, 41661u, 0 },
   { 2746102297u, 162087u, 0 },
   { 2746102297u, 162087u, 0 },
   { 1936779721u, 664710u, 0 },
   { 1070451441u, 1501768u, 0 },
   { 1070451441u, 1501768u, 0 },
   { 2061289617u, 2710474u, 0 },
   { 2061289617u, 2710474u, 0 },
   { 4235760785u, 44382509u, 0 },
   { 2312776601u, 45783875u, 0 },
   { 2678348049u, 165920782u, 0 },
   { 3315991761u, 413007985u, 0 },
   { 1567179849u, 541956877u, 0 },
   { 2273104657u, 1486621767u, 0 },
   { 3796117489u, 1867116582u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 3763487577u, 3377920039u, 3u },
   { 2972093785u, 1402148275u, 11u },
   { 2785759393u, 3968325740u, 28u },
   { 551239041u, 3335735663u, 50u },
   { 551239041u, 3335735663u, 50u },
   { 732515297u, 554264481u, 116u },
   { 732515297u, 554264481u, 116u },
   { 732515297u, 554264481u, 116u },
   { 2681625681u, 3960593856u, 739u },
   { 2546105329u, 1679561907u, 1875u },
   { 533319201u, 2248012685u, 5393u },
   { 533319201u, 2248012685u, 5393u },
   { 996692113u, 2949507147u, 16011u },
   { 996692113u, 2949507147u, 16011u },
   { 2616068761u, 328479117u, 198156u },
   { 1411295841u, 761797252u, 229581u } 
};
#else
mp_limb_t flint_fmpz_pseudosquares[][2] = 
{
   { 17, 0 },
   { 73, 0 },
   { 241, 0 },
   { 1009, 0 },
   { 2641, 0 },
   { 8089, 0 },
   { 18001, 0 },
   { 53881, 0 },
   { 87481, 0 },
   { 117049, 0 },
   { 515761, 0 },
   { 1083289, 0 },
   { 3206641, 0 },
   { 3818929, 0 },
   { 9257329, 0 },
   { 22000801, 0 },
   { 48473881, 0 },
   { 48473881, 0 },
   { 175244281, 0 },
   { 427733329, 0 },
   { 427733329, 0 },
   { 898716289u, 0 },
   { 2805544681u, 0 },
   { 2805544681u, 0 },
   { 2805544681u, 0 },
   { 10310263441u, 0 },
   { 23616331489u, 0 },
   { 85157610409u, 0 },
   { 85157610409u, 0 },
   { 196265095009u, 0 },
   { 196265095009u, 0 },
   { 2871842842801u, 0 },
   { 2871842842801u, 0 },
   { 2871842842801u, 0 },
   { 26250887023729u, 0 },
   { 26250887023729u, 0 },
   { 112434732901969u, 0 },
   { 112434732901969u, 0 },
   { 112434732901969u, 0 },
   { 178936222537081u, 0 },
   { 178936222537081u, 0 },
   { 696161110209049u, 0 },
   { 696161110209049u, 0 },
   { 2854909648103881u, 0 },
   { 6450045516630769u, 0 },
   { 6450045516630769u, 0 },
   { 11641399247947921u, 0 },
   { 11641399247947921u, 0 },
   { 190621428905186449u, 0 },
   { 196640248121928601u, 0 },
   { 712624335095093521u, 0 },
   { 1773855791877850321u, 0 },
   { 2327687064124474441u, 0 },
   { 6384991873059836689u, 0 },
   { 8019204661305419761u, 0 },
   { 10198100582046287689u, 0 },
   { 10198100582046287689u, 0 },
   { 10198100582046287689u, 0 },
   { 14508056099771532121u, 3u },
   { 6022180988239908185u, 11u },
   { 17043829275960758433u, 28u },
   { 14326875581237116289u, 50u },
   { 14326875581237116289u, 50u },
   { 2380547819961928673u, 116u },
   { 2380547819961928673u, 116u },
   { 2380547819961928673u, 116u },
   { 17010621086940159057u, 739u },
   { 7213663464718498801u, 1875u },
   { 9655140963601468961u, 5393u },
   { 9655140963601468961u, 5393u },
   { 12668036736679956625u, 16011u },
   { 12668036736679956625u, 16011u },
   { 1410807067550026393u, 198156u },
   { 3271894284933966433u, 229581u } 
};
#endif

#define FLINT_NUM_FMPZ_PSEUDOSQUARES 74

void fmpz_set_pseudosquare(fmpz_t f, unsigned int i)
{
#ifndef FLINT64
   if (i < 25)
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][0]);
   else if (i < 58)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, flint_fmpz_pseudosquares[i][0]);
   } else if (i < FLINT_NUM_FMPZ_PSEUDOSQUARES)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][2]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, flint_fmpz_pseudosquares[i][1]);
   } 
#else
   if (i < 58)
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][0]);
   else if (i < FLINT_NUM_FMPZ_PSEUDOSQUARES)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 64);
      fmpz_add_ui(f, f, flint_fmpz_pseudosquares[i][0]);
   }      
#endif
   else
   {
      printf("Exception (fmpz_set_pseudosquare). Index too large.\n");
      abort();
   }
}

int fmpz_is_prime_pseudosquare(fmpz_t n)
{
    unsigned int i, j, m1;
    mp_limb_t p, B, mod8;
    fmpz_t NB, f, exp, mod, nm1;
    int ret;

    if (fmpz_sgn(n) <= 0) 
       return 0;

    if (fmpz_size(n) == 1) 
       return n_is_prime_pseudosquare(fmpz_get_ui(n));

    n_compute_primes(FLINT_PSEUDOSQUARES_CUTOFF);

    for (i = 0; i < FLINT_PSEUDOSQUARES_CUTOFF; i++)
    {
        p = flint_primes[i];
        if (fmpz_fdiv_ui(n, p) == 0) 
           return 0;
    }

    fmpz_init(NB);
    fmpz_init(f);
    fmpz_init(exp);
    fmpz_init(mod);
    fmpz_init(nm1);
    
    B  = flint_primes[FLINT_PSEUDOSQUARES_CUTOFF];
    fmpz_sub_ui(nm1, n, 1);
    fmpz_fdiv_q_ui(NB, nm1, B);
    fmpz_add_ui(NB, NB, 1);
    
    m1 = 0;

    for (i = 0; i < FLINT_NUM_FMPZ_PSEUDOSQUARES; i++)
    {
       fmpz_set_pseudosquare(f, i);
       if (fmpz_cmp(f, NB) > 0) 
          break;
    }

    if (i == FLINT_NUM_FMPZ_PSEUDOSQUARES)
    {
       ret = -1;
       goto cleanup;
    }

    fmpz_fdiv_q_2exp(exp, nm1, 1);
    
    for (j = 0; j <= i; j++)
    {
        fmpz_set_ui(mod, flint_primes[j]);
        fmpz_powm(mod, mod, exp, n);
        if (!fmpz_is_one(mod) && fmpz_cmp(mod, nm1) != 0) 
        {
           ret = 0;
           goto cleanup;
        }
        if (fmpz_cmp(mod, nm1) == 0) 
           m1 = 1;
    }

    mod8 = fmpz_fdiv_ui(n, 8);

    if ((mod8 == 3) || (mod8 == 7)) 
    {
       ret = 1;
       goto cleanup;
    }
            
    if (mod8 == 5)
    {
        fmpz_set_ui(mod, 2);
        fmpz_powm(mod, mod, exp, n);
        if (fmpz_cmp(mod, nm1) == 0)
        {
           ret = 1;
           goto cleanup;
        }
        printf("Whoah, ");
        fmpz_print(n);
        printf("is a probable prime, but not prime, please report!!\n");
        abort();
    }
    else
    {
        if (m1) 
        {
            ret = 1;
            goto cleanup;
        }
            
        for (j = i + 1; j < FLINT_NUM_FMPZ_PSEUDOSQUARES + 1; j++)
        {
            fmpz_set_ui(mod, flint_primes[j]);
            fmpz_powm(mod, mod, exp, n);
            if (fmpz_cmp(mod, nm1) == 0)
            {
               ret = 1;
               goto cleanup;
            }
            if (!fmpz_is_one(mod))
            {
                printf("Whoah, ");
                fmpz_print(n);
                printf("is a probable prime, but not prime, please report!!\n");
                abort();
            }
        }
        printf("Whoah, ");
        fmpz_print(n);
        printf("is a probable prime, but not prime, please report!!\n");
        abort();
    }

cleanup:

    fmpz_clear(NB);
    fmpz_clear(f);
    fmpz_clear(exp);
    fmpz_clear(mod);
    fmpz_clear(nm1);

    return ret;
}
