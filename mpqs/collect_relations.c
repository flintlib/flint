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

   Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#include <gmp.h>
#include <string.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

void
mpqs_do_sieving(mpqs_t mpqs_inf, unsigned char * sieve)
{
    slong num_primes = mpqs_inf->num_primes;
    mp_limb_t * soln1 = mpqs_inf->soln1;
    mp_limb_t * soln2 = mpqs_inf->soln2;
    prime_t * factor_base = mpqs_inf->factor_base;
    mp_limb_t p;
    unsigned char * end = sieve + mpqs_inf->sieve_size;
    register unsigned char * pos1;
    register unsigned char * pos2;
    register unsigned char * bound;
    slong size;
    slong diff;
    slong pind;

    memset(sieve, 0, mpqs_inf->sieve_size + sizeof(ulong));
    *end = (char) 255;

    for (pind = mpqs_inf->small_primes; pind < num_primes; pind++)
    {
        p = factor_base[pind].p;
        size = factor_base[pind].size;
        pos1 = sieve + soln1[pind];
        pos2 = sieve + soln2[pind];
        diff = pos2 - pos1;
        bound = end - 2*p;

        while (bound - pos1 > 0)
        {
            (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
            (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
        }

        while ((end - pos1 > 0) && (end - pos1 - diff > 0))
            (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;

        pos2 = pos1 + diff;

        if (end - pos2 > 0)
            (*pos2) += size;

        if (end - pos1 > 0)
            (*pos1) += size;
    }
}

slong
mpqs_evaluate_candidate(mpqs_t mpqs_inf, slong i, unsigned char * sieve)
{
    slong bits, exp, extra_bits;
    mp_limb_t modp, prime;
    slong num_primes = mpqs_inf->num_primes;
    prime_t * factor_base = mpqs_inf->factor_base;
    fac_t * factor = mpqs_inf->factor;
    mp_limb_t * soln1 = mpqs_inf->soln1;
    mp_limb_t * soln2 = mpqs_inf->soln2;
    slong * small = mpqs_inf->small;
    mp_limb_t pinv;
    slong num_factors = 0;
    slong relations = 0;
    slong j;

    fmpz_t X, Y, res, p;
    fmpz_init(X);
    fmpz_init(Y);
    fmpz_init(res);
    fmpz_init(p);

    fmpz_set_ui(X, i - mpqs_inf->sieve_size / 2); /* X */

    fmpz_mul(Y, X, mpqs_inf->A);
    fmpz_add(Y, Y, mpqs_inf->B); /* Y = AX+B */
    fmpz_add(res, Y, mpqs_inf->B);

    fmpz_mul(res, res, X);
    fmpz_add(res, res, mpqs_inf->C); /* res = AX^2 + 2BX + C */

    bits = FLINT_ABS(fmpz_bits(res));
    bits -= BITS_ADJUST;
    extra_bits = 0;

    if (factor_base[0].p != 1) /* divide out powers of the multiplier */
    {
      fmpz_set_ui(p, factor_base[0].p);
      exp = fmpz_remove(res, res, p);
      if (exp) extra_bits += exp*mpqs_inf->factor_base[0].size;
      small[0] = exp;
    }
    else small[0] = 0;

    fmpz_set_ui(p, 2); /* divide out by powers of 2 */
    exp = fmpz_remove(res, res, p);

    extra_bits += exp;
    small[1] = exp;

    for (j = 2; j < mpqs_inf->small_primes; j++) /* pull out small primes */
    {
        prime = factor_base[j].p;
        pinv = factor_base[j].pinv;
        modp = n_mod2_preinv(i, prime, pinv);

        if ((modp == soln1[j]) || (modp == soln2[j]))
        {
            fmpz_set_ui(p, prime);
            exp = fmpz_remove(res, res, p);
            if (exp) extra_bits += mpqs_inf->factor_base[j].size;
            small[j] = exp;
        }
        else small[j] = 0;
    }

    if (extra_bits + sieve[i] > bits)
    {
        sieve[i] += extra_bits;

        /* pull out remaining primes */
        for (j = mpqs_inf->small_primes; j < num_primes && extra_bits < sieve[i]; j++)
        {
            prime = factor_base[j].p;
            pinv = factor_base[j].pinv;
            modp = n_mod2_preinv(i, prime, pinv);

            if (soln2[j] != 0)
            {
                if ((modp == soln1[j]) || (modp == soln2[j]))
                {
                    fmpz_set_ui(p, prime);
                    exp = fmpz_remove(res, res, p);
                    if (exp)
                    {
                        extra_bits += mpqs_inf->factor_base[j].size;
                        factor[num_factors].ind = j;
                        factor[num_factors++].exp = exp;
                    }
                }
            }
            else
            {
                fmpz_set_ui(p, prime);
                exp = fmpz_remove(res, res, p);
                factor[num_factors].ind = j;
                factor[num_factors++].exp = exp + 1;
            }
        }

        if (fmpz_cmp_ui(res, 1) == 0 || fmpz_cmp_si(res, -1) == 0) /* We've found a relation */
        {            
            mpqs_inf->num_factors = num_factors;

            relations += mpqs_insert_relation(mpqs_inf, Y);  /* Insert the relation in the matrix */

            if (mpqs_inf->num_relations >= mpqs_inf->buffer_size)
            {
                flint_printf("Error: too many duplicate relations!\n");
                abort();
            }

            goto cleanup;
        }
    }

    cleanup:
    fmpz_clear(X);
    fmpz_clear(Y);
    fmpz_clear(res);
    fmpz_clear(p);

    return relations;
}


slong
mpqs_evaluate_sieve(mpqs_t mpqs_inf, unsigned char * sieve)
{
    slong i = 0, j = 0;
    ulong * sieve2 = (ulong *) sieve;
    unsigned char bits = mpqs_inf->sieve_bits;
    slong rels = 0;

    while (j < mpqs_inf->sieve_size / sizeof(ulong))
    {

#if FLINT64
        while ((sieve2[j] & UWORD(0xC0C0C0C0C0C0C0C0)) == 0)
#else
        while ((sieve2[j] & UWORD(0xC0C0C0C0)) == 0)
#endif
        {
            j++;
        }

        i = j * sizeof(ulong);

        while (i < (j + 1) * sizeof(ulong) && i < mpqs_inf->sieve_size)
        {
            if (sieve[i] > bits)
                rels += mpqs_evaluate_candidate(mpqs_inf, i, sieve);

            i++;
        }
        j++;
    }

    return rels;
}

slong
mpqs_collect_relations(mpqs_t mpqs_inf, unsigned char * sieve)
{
    slong relations = 0;

    mpqs_do_sieving(mpqs_inf, sieve);
    relations += mpqs_evaluate_sieve(mpqs_inf, sieve);
   
    return relations;
}
