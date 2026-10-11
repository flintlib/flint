/*
    Copyright (C) 2006, 2011, 2016, 2020 William Hart
    Copyright (C) 2015 Nitin Kumar
    Copyright (C) 2020 Dan Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>
#include "thread_pool.h"
#include <gmp.h>
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

/*
    The actual sieving part of the quadratic sieve. This version is only run
    if the sieve block is small, i.e. less than 2*BLOCK_SIZE.
*/
void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
   slong num_primes = qs_inf->num_primes;
   int * soln1 = poly->soln1;
   int * soln2 = poly->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   ulong p;

   unsigned char * end = sieve + qs_inf->sieve_size;
   register unsigned char * pos1;
   register unsigned char * pos2;
   register unsigned char * bound;
   slong size;
   slong diff;
   slong pind;

   /* set entries in sieve to initial value and put a sentinel at the end */
   memset(sieve, qs_inf->sieve_fill, qs_inf->sieve_size + sizeof(ulong));
   *end = (char) 255;

   for (pind = qs_inf->small_primes; pind < num_primes; pind++)
   {
      if (soln2[pind] == 0) continue; /* don't sieve with A factors */

      p = factor_base[pind].p;
      size = factor_base[pind].size;
      pos1 = sieve + soln1[pind];
      pos2 = sieve + soln2[pind];
      diff = pos2 - pos1;
      bound = end - 2*p;

      /* mark off entries in sieve, unrolled by 2 */
      while (bound - pos1 > 0)
      {
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
      }

      /* deal with final entries at the end, missed by unrolled loop */
      while ((end - pos1 > 0) && (end - pos1 - diff > 0))
      {
         (*pos1) += size, (*(pos1 + diff)) += size, pos1 += p;
      }

      pos2 = pos1 + diff;

      if (end - pos2 > 0)
      {
         (*pos2) += size;
      }

      if (end - pos1 > 0)
      {
         (*pos1) += size;
      }
   }
}

/*
   Sieve the block [., B) with the factor base primes of index start to
   stop - 1, continuing from (and updating) the positions posn1/posn2.
   posn1[i] is the next sieve position hit by prime i, and posn2[i] the
   distance from there to the next hit by the other root.  Nothing assumes
   the primes are small compared with the block, but the unrolled loop is
   only effective for primes that hit the block several times.
*/
static inline void
qsieve_sieve_block(unsigned char * sieve, unsigned char * B, slong start,
     slong stop, const prime_t * factor_base, const int * soln2,
     int * posn1, int * posn2)
{
    slong pind, d1, d2, size;
    ulong p;
    unsigned char * Bp;
    unsigned char * pos;

    for (pind = start; pind < stop; pind++)
    {
        if (soln2[pind] == 0) /* skip primes dividing A */
            continue;

        p = factor_base[pind].p;
        size = factor_base[pind].size;
        d1 = posn2[pind];
        d2 = p - d1;
        Bp = B - 2*d1 - d2;
        pos = sieve + posn1[pind];

        while (pos < Bp)
        {
            (*pos) += size, (*(pos + d1)) += size, pos += p;
            (*pos) += size, (*(pos + d1)) += size, pos += p;
        }

        Bp = B - d1;

        while (pos < Bp)
        {
            (*pos) += size,
            (*(pos + d1)) += size, pos += p;
        }

        if (pos < B)
        {
            (*pos) += size, pos += d1;
            posn2[pind] = d2;
        }
        else
        {
            posn2[pind] = d1;
        }

        posn1[pind] = (pos - sieve);
    }
}

/*
   Second sieving routine. This version breaks sieve array into blocks of
   size BLOCK_SIZE, and sieves each block with the whole factor base.  Within
   each block, the primes below QS_SMALL_BOUND are sieved over sub-blocks of
   QS_L1_SIZE bytes, which fit in L1 cache; these primes account for most of
   the sieve updates.  A final partial block is handled, so the sieve size
   need not be a multiple of BLOCK_SIZE (it used to be left unsieved).
*/
void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong i, start, end, sub, subend, pind, l1_prime, second_prime, size;
    ulong p;
    slong num_primes = qs_inf->num_primes;
    slong sieve_size = qs_inf->sieve_size;
    int * soln1 = poly->soln1;
    int * soln2 = poly->soln2;
    int * posn1 = poly->posn1;
    int * posn2 = poly->posn2;
    prime_t * factor_base = qs_inf->factor_base;

    unsigned char * B;
    register unsigned char * pos;

    /* fill sieve with initial values and put sentinel at end */
    memset(sieve, qs_inf->sieve_fill, sieve_size + sizeof(ulong));
    sieve[sieve_size] = (char) 255;

    /*
       initial values for positions (which must be saved at the end of each
       sieve block in preparation for start of next sieve block)
    */
    for (i = 0; i < num_primes; i++)
    {
        posn1[i] = soln1[i];
        posn2[i] = soln2[i] - posn1[i];
    }

    second_prime = qs_inf->second_prime;

    /* index of the first prime that is not sieved in L1 sized blocks */
    for (l1_prime = qs_inf->small_primes; l1_prime < second_prime; l1_prime++)
    {
        if (factor_base[l1_prime].p >= QS_SMALL_BOUND)
            break;
    }

    for (start = 0; start < sieve_size; start += BLOCK_SIZE)
    {
        /* end of current sieve block */
        end = FLINT_MIN(start + BLOCK_SIZE, sieve_size);
        B = sieve + end;

        /* small primes, which hit an L1 sized block many times */
        for (sub = start; sub < end; sub += QS_L1_SIZE)
        {
            subend = FLINT_MIN(sub + QS_L1_SIZE, end);
            qsieve_sieve_block(sieve, sieve + subend, qs_inf->small_primes,
                         l1_prime, factor_base, soln2, posn1, posn2);
        }

        /* medium sized primes, which hit the sieve block at least once */
        qsieve_sieve_block(sieve, B, l1_prime, second_prime,
                         factor_base, soln2, posn1, posn2);

        /*
            now deal with larger primes which are likely to only hit sieve
            block once, if at all
        */
        for (pind = second_prime; pind < num_primes; pind++)
        {
            p = factor_base[pind].p;

            if (soln2[pind] == 0)
                continue;

            size = factor_base[pind].size;
            pos = sieve + posn1[pind];

            if (pos < B) /* hits the sieve block */
            {
                (*pos) += size;
                pos += posn2[pind];

                if (pos < B)
                {
                    (*pos) += size;
                    pos += p - posn2[pind];
                } else
                {
                    posn2[pind] = p - posn2[pind];
                }

                posn1[pind] = (pos - sieve);
            }
        }
    }
}


#define QS_BUCKET_CHUNK 256 /* primes added to the buckets between checks */

/* make room for at least `need` more entries in bucket b */
static void
qsieve_bucket_grow(qs_poly_t poly, slong b, slong need)
{
    slong used = poly->bucket_ptr[b] - poly->bucket[b];
    slong alloc = FLINT_MAX(2*poly->bucket_alloc[b], used + need);

    poly->bucket[b] = flint_realloc(poly->bucket[b], alloc*sizeof(uint32_t));
    poly->bucket_ptr[b] = poly->bucket[b] + used;
    poly->bucket_alloc[b] = alloc;
}

/*
   Third sieving routine, as the second one but with a bucket sieve for the
   primes larger than BLOCK_SIZE; see the comment at BLOCK_SIZE in qsieve.h.
*/
void qsieve_do_sieving3(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong b, big, i, j, k, j0, j1, jj, jend, start, stop, m;
    slong l1_prime = qs_inf->l1_prime;
    slong num_blocks = qs_inf->num_blocks;
    slong num_slices = qs_inf->num_slices;
    slong bucket_prime = qs_inf->bucket_prime;
    slong sieve_size = qs_inf->sieve_size;
    slong * slice_start = qs_inf->slice_start;
    slong * slice_end = poly->slice_end;
    int * soln1 = poly->soln1;
    int * soln2 = poly->soln2;
    int * posn1 = poly->posn1;
    int * posn2 = poly->posn2;
    uint32_t ** bucket = poly->bucket;
    uint32_t ** bucket_ptr = poly->bucket_ptr;
    prime_t * factor_base = qs_inf->factor_base;
    const uint32_t * p32 = qs_inf->fb_p32;
    uint32_t size32 = sieve_size;
    ulong p;
    unsigned char size;

    /* fill sieve with initial values and put sentinel at end */
    memset(sieve, qs_inf->sieve_fill, sieve_size + sizeof(ulong));
    sieve[sieve_size] = (char) 255;

    /* fill the buckets */
    for (b = 0; b < num_blocks + 2; b++)
        bucket_ptr[b] = bucket[b];

    for (k = 0; k < num_slices; k++)
    {
        j0 = slice_start[k];
        j1 = slice_start[k + 1];

        for (j = j0; j < j1; j = jend)
        {
            slong need;

            jend = FLINT_MIN(j + QS_BUCKET_CHUNK, j1);

            /* bound on the hits per block of the primes in this chunk */
            p = factor_base[j].p;
            need = 2*(jend - j)*((QS_L1_SIZE + p - 1) / p);

            for (b = 0; b < num_blocks; b++)
            {
                if (poly->bucket_alloc[b] - (bucket_ptr[b] - bucket[b]) < need)
                    qsieve_bucket_grow(poly, b, need);
            }

            /*
               If each root hits the sieve at most m <= 4 times, where m is
               the same for the whole chunk, do all m possible hits without
               branching on whether they are in the sieve: hits outside it
               are written to one of two dummy buckets, which are never read.
               The loops then have a fixed number of iterations, so they
               avoid the mispredicted branches which otherwise dominate the
               cost for these primes.  (The primes dividing A are much
               smaller than these.)
            */
            m = (sieve_size + p - 1) / p;

            if (j >= qs_inf->huge_prime || m <= 4)
            {
                if (j >= qs_inf->huge_prime)
                    m = 1;

                for (jj = j; jj < jend; jj++)
                {
                    uint32_t hi, r1, r2, b1, b2, pp, t;
                    uint32_t * q;

                    r1 = soln1[jj];
                    r2 = soln2[jj];
                    pp = p32[jj];
                    hi = ((uint32_t) (jj - j0)) << 16;

                    for (t = 0; t < m; t++)
                    {
                        b1 = (r1 < size32) ? (r1 >> QS_L1_BITS) : num_blocks;
                        b2 = (r2 < size32) ? (r2 >> QS_L1_BITS) : num_blocks + 1;

                        q = bucket_ptr[b1];
                        *q = hi | (r1 & (QS_L1_SIZE - 1));
                        bucket_ptr[b1] = q + 1;
                        q = bucket_ptr[b2];
                        *q = hi | (r2 & (QS_L1_SIZE - 1));
                        bucket_ptr[b2] = q + 1;

                        r1 += pp;
                        r2 += pp;
                    }
                }

                bucket_ptr[num_blocks] = bucket[num_blocks];
                bucket_ptr[num_blocks + 1] = bucket[num_blocks + 1];
                continue;
            }

            for (jj = j; jj < jend; jj++)
            {
                uint32_t hi, r1, r2, pp;
                uint32_t * q;

                r1 = soln1[jj];
                r2 = soln2[jj];

                if (r2 == 0) /* skip primes dividing A */
                    continue;

                pp = p32[jj];
                hi = ((uint32_t) (jj - j0)) << 16;

                /* r1 < r2, so r1 is in range whenever r2 is */
                while (r2 < size32)
                {
                    q = bucket_ptr[r1 >> QS_L1_BITS];
                    *q = hi | (r1 & (QS_L1_SIZE - 1));
                    bucket_ptr[r1 >> QS_L1_BITS] = q + 1;
                    q = bucket_ptr[r2 >> QS_L1_BITS];
                    *q = hi | (r2 & (QS_L1_SIZE - 1));
                    bucket_ptr[r2 >> QS_L1_BITS] = q + 1;
                    r1 += pp;
                    r2 += pp;
                }

                if (r1 < size32)
                {
                    q = bucket_ptr[r1 >> QS_L1_BITS];
                    *q = hi | (r1 & (QS_L1_SIZE - 1));
                    bucket_ptr[r1 >> QS_L1_BITS] = q + 1;
                }
            }
        }

        for (b = 0; b < num_blocks; b++)
            slice_end[k*num_blocks + b] = bucket_ptr[b] - bucket[b];
    }

    /* positions of the primes sieved directly */
    for (i = qs_inf->small_primes; i < bucket_prime; i++)
    {
        posn1[i] = soln1[i];
        posn2[i] = soln2[i] - posn1[i];
    }

    for (big = 0; big < sieve_size; big += BLOCK_SIZE)
    {
        slong bigend = FLINT_MIN(big + BLOCK_SIZE, sieve_size);

        for (b = big / QS_L1_SIZE; b*QS_L1_SIZE < bigend; b++)
        {
            unsigned char * block = sieve + b*QS_L1_SIZE;
            slong end = FLINT_MIN((b + 1)*QS_L1_SIZE, sieve_size);

            /* small primes, directly over the L1 sized block */
            qsieve_sieve_block(sieve, sieve + end, qs_inf->small_primes,
                            l1_prime, factor_base, soln2, posn1, posn2);

            /* large primes, from the bucket of the block */
            start = 0;
            for (k = 0; k < num_slices; k++)
            {
                uint32_t * e = bucket[b];

                stop = slice_end[k*num_blocks + b];
                size = qs_inf->slice_size[k];

                for (i = start; i < stop; i++)
                    block[e[i] & 0xffff] += size;

                start = stop;
            }
        }

        /* medium primes, directly over the whole block */
        qsieve_sieve_block(sieve, sieve + bigend, l1_prime, bucket_prime,
                            factor_base, soln2, posn1, posn2);
    }
}

/*
    Remove all powers of p from the two-word value {hi,lo}, returning the
    exponent.  Sieve residues fit in two words for inputs up to roughly 200
    bits, where this is about twice as fast as fmpz_remove.

    For odd p the quotient comes from exact division: p divides exactly, so
    one multiplication by p^(-1) mod 2^FLINT_BITS per limb suffices.  p = 2 is
    a shift.  The multiplier held in factor_base[0] may be even and composite,
    which neither case covers, so that falls back to mpn.
*/
static ulong
qsieve_ll_remove(ulong * hip, ulong * lop, ulong p, ulong pinv)
{
    ulong exp = 0, hi = *hip, lo = *lop;

    if ((hi | lo) == 0)
        return 0;

    if (p == 2)
    {
        if (lo == 0)
        {
            exp = FLINT_BITS + flint_ctz(hi);
            lo = hi >> flint_ctz(hi);
            hi = 0;
        }
        else
        {
            exp = flint_ctz(lo);
            if (exp)
            {
                lo = (lo >> exp) | (hi << (FLINT_BITS - exp));
                hi >>= exp;
            }
        }
    }
    else if (p & 1)
    {
        ulong pbinv = n_binvert(p);

        while (1)
        {
            if (hi == 0)
            {
                ulong q = lo / p;

                if (q * p != lo)
                    break;

                lo = q;
                exp++;

                if (lo == 1)
                    break;
            }
            else
            {
                ulong q0, c;

                if (n_ll_mod_preinv(hi, lo, p, pinv) != 0)
                    break;

                q0 = lo * pbinv;
                c = n_mulhi(q0, p);
                lo = q0;
                hi = (hi - c) * pbinv;
                exp++;
            }
        }
    }
    else
    {
        ulong t[2];

        while (1)
        {
            slong n = hi ? 2 : 1;

            t[0] = lo;
            t[1] = hi;

            if (mpn_mod_1(t, n, p) != 0)
                break;

            mpn_divexact_1(t, t, n, p);
            lo = t[0];
            hi = (n == 2) ? t[1] : 0;
            exp++;

            if (hi == 0 && lo == 1)
                break;
        }
    }

    *hip = hi;
    *lop = lo;
    return exp;
}

/* sort factors by index (insertion sort, the lists are short) */
static void
qsieve_sort_factors(fac_t * factor, slong num)
{
    slong i, j;
    fac_t t;

    for (i = 1; i < num; i++)
    {
        t = factor[i];
        for (j = i; j > 0 && factor[j - 1].ind > t.ind; j--)
            factor[j] = factor[j - 1];
        factor[j] = t;
    }
}

/* remove p from the residue, whichever representation is in use */
#define QS_REMOVE(prime_, preinv_)                                          \
    (use_ll ? qsieve_ll_remove(&rhi, &rlo, (prime_), (preinv_))             \
            : (fmpz_set_ui(p, (prime_)), fmpz_remove(res, res, p)))

/*
    check position i in sieve array for smoothness
*/
slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)
{
   slong bits, exp, extra_bits;
   ulong modp, prime;
   slong num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   slong * small = poly->small; /* exponents of small primes and mult. */
   fac_t * factor = poly->factor;
   int * soln1 = poly->soln1;
   int * soln2 = poly->soln2;
   ulong * A_ind = qs_inf->A_ind;
   ulong pinv;
   slong num_factors = 0;
   slong relations = 0;
   slong j, k, top;

   fmpz_t X, Y, res, p;
   ulong rhi = 0, rlo = 0;
   int use_ll, neg = 0;

   fmpz_init(X);
   fmpz_init(Y);
   fmpz_init(res);
   fmpz_init(p);

   /* C = (B^2 - kn)/A only depends on the polynomial */
   if (!poly->C_valid)
   {
      qsieve_compute_C(poly->C, qs_inf, poly);

      /*
         Modulo a prime q dividing A, the polynomial A*x^2 + 2*B*x + C is
         linear, as 2B is invertible mod q (B^2 = kn mod q), so it has the
         single root x = -C/(2B).  Record it as a sieve position, so that
         the candidates can be checked for q cheaply.
      */
      for (k = 0; k < qs_inf->s; k++)
      {
         ulong q = factor_base[A_ind[k]].p, b2, cq, x;

         b2 = n_mulmod2_preinv(fmpz_fdiv_ui(poly->B, q), 2, q, factor_base[A_ind[k]].pinv);
         cq = fmpz_fdiv_ui(poly->C, q);
         x = n_mulmod2_preinv(n_negmod(cq, q), n_invmod(b2, q), q, factor_base[A_ind[k]].pinv);
         poly->A_root[k] = n_addmod(x, (qs_inf->sieve_size / 2) % q, q);
      }

      poly->C_valid = 1;
   }

   fmpz_set_si(X, i - qs_inf->sieve_size / 2); /* X */

   fmpz_mul(Y, X, qs_inf->A);
   fmpz_add(Y, Y, poly->B); /* Y = AX+B */
   fmpz_add(res, Y, poly->B); /* Y = AX+2B */

   fmpz_mul(res, res, X);
   fmpz_add(res, res, poly->C); /* res = AX^2 + 2BX + C */

#if QS_DEBUG & 128
   flint_printf("res = "); fmpz_print(res); flint_printf("\n");
   flint_printf("Poly: "); fmpz_print(qs_inf->A); flint_printf("*x^2 + 2*");
   fmpz_print(poly->B); flint_printf("*x + "); fmpz_print(C); flint_printf("\n");
   flint_printf("x = %wd\n", i - qs_inf->sieve_size / 2);
#endif

   sieve[i] -= qs_inf->sieve_fill; /* adjust sieve entry to number of bits */
   bits = FLINT_ABS(fmpz_bits(res)); /* compute bits of poly value */
   bits -= BITS_ADJUST; /* adjust for log approximations */
   extra_bits = 0; /* bits for mult. and small primes we didn't sieve with */

   /*
      Work in two words when the polynomial value fits, which it does for all
      but the largest inputs.  The sign is taken out here and put back in
      small[2] below, exactly as the fmpz path does at the end.
   */
   /*
      Only worth it when the value really needs two words.  Below about 100
      bits the polynomial values fit in one, where fmpz_remove already stays
      on its own single limb fast path and this would just add a branch.
   */
   use_ll = 0;

   if (COEFF_IS_MPZ(*res))
   {
      mpz_ptr z = COEFF_TO_PTR(*res);
      slong sz = z->_mp_size;

      if (FLINT_ABS(sz) == 2)
      {
         use_ll = 1;
         neg = (sz < 0);
         rlo = z->_mp_d[0];
         rhi = z->_mp_d[1];
      }
   }

   if (factor_base[0].p != 1) /* divide out powers of the multiplier */
   {
      exp = QS_REMOVE(factor_base[0].p, factor_base[0].pinv);
      if (exp)
          extra_bits += exp*qs_inf->factor_base[0].size;
      small[0] = exp;
#if QS_DEBUG & 128
      if (exp)
          flint_printf("%ld^%ld ", factor_base[0].p, exp);
#endif
   }
   else
   {
       small[0] = 0;
   }

   exp = QS_REMOVE(UWORD(2), factor_base[1].pinv); /* divide out powers of 2 */
#if QS_DEBUG & 128
   if (exp)
       flint_printf("%ld^%ld ", 2, exp);
#endif

   extra_bits += exp;
   small[1] = exp;

   for (j = 3; j < qs_inf->small_primes; j++) /* pull out small primes */
   {
      prime = factor_base[j].p;
      pinv = factor_base[j].pinv2;
      FLINT_ASSERT(prime < UWORD(1) << (FLINT_BITS / 2));
      modp = n_mod_lemire(i, prime, pinv);

      if (modp == soln1[j] || modp == soln2[j])
      {
         exp = QS_REMOVE(prime, factor_base[j].pinv);
         if (exp)
             extra_bits += qs_inf->factor_base[j].size;
         small[j] = exp;
#if QS_DEBUG & 128
         if (exp)
            flint_printf("%ld^%ld ", prime, exp);
#endif
      }
      else
      {
         small[j] = 0;
      }
   }

   /* if we have a chance of being a candidate */
   if (extra_bits + sieve[i] > bits)
   {
      sieve[i] += extra_bits; /* add bits from small primes and mult. */

      /*
         The primes dividing A are skipped by the sieve, so they do not
         count towards the sieve value.  Remove them first, so that the
         search below can stop as soon as the primes found account for the
         sieve value without missing any of them.
      */
      for (k = 0; k < qs_inf->s; k++)
      {
         j = A_ind[k];

         if (n_mod2_preinv(i, factor_base[j].p, factor_base[j].pinv) == poly->A_root[k])
            exp = QS_REMOVE(factor_base[j].p, factor_base[j].pinv);
         else
            exp = 0;

         factor[num_factors].ind = j;
         factor[num_factors++].exp = exp + 1; /* really factoring A*f(i) */
      }

      top = num_primes;

      if (qs_inf->use_buckets)
      {
         /*
            The bucket primes dividing the polynomial value are exactly those
            with an entry for position i in the bucket of its block, so read
            them off instead of trial dividing.
         */
         slong blk = i >> QS_L1_BITS, start = 0, stop, t;
         uint32_t off = i & (QS_L1_SIZE - 1);
         uint32_t * e = poly->bucket[blk];

         top = qs_inf->bucket_prime;

         for (k = 0; k < qs_inf->num_slices; k++)
         {
            stop = poly->slice_end[k*qs_inf->num_blocks + blk];

            /* look for matches a chunk at a time, which vectorises */
            for (t = start; t < stop; t += QS_TRIAL_CHUNK)
            {
               slong tend = FLINT_MIN(t + QS_TRIAL_CHUNK, stop), u;
               unsigned char hit[QS_TRIAL_CHUNK];
               int any = 0;

               for (u = t; u < tend; u++)
               {
                  hit[u - t] = ((uint16_t) e[u] == (uint16_t) off);
                  any |= hit[u - t];
               }

               if (!any)
                  continue;

               for (u = t; u < tend; u++)
               {
                  if (hit[u - t])
                  {
                     j = qs_inf->slice_start[k] + (e[u] >> 16);

                     if (num_factors >= qs_inf->max_factors)
                        goto cleanup;

                     exp = QS_REMOVE(factor_base[j].p, factor_base[j].pinv);
                     extra_bits += factor_base[j].size;
                     factor[num_factors].ind = j;
                     factor[num_factors++].exp = exp;
                  }
               }
            }

            start = stop;
         }
      }

      /*
         Trial divide by the remaining primes, a chunk at a time, stopping
         once the primes found account for the sieve value.  Prime j
         divides the value at position i if and only if i - soln1[j] or
         i - soln2[j] is divisible by p_j.  For odd p < 2^32 and 0 <= x <
         2^32, p divides x if and only if x * p^(-1) mod 2^32 <= (2^32-1)/p,
         so this can be checked with word size multiplications, which the
         compiler can vectorise.  Primes dividing A have soln1[j] = soln2[j]
         = 0 and were dealt with above, so they are skipped.
      */
      {
         const uint32_t * p32 = qs_inf->fb_p32;
         const uint32_t * inv32 = qs_inf->fb_inv32;
         const uint32_t * lim32 = qs_inf->fb_lim32;
         unsigned char hit[QS_TRIAL_CHUNK];
         slong jend, t;
         uint32_t i32 = (uint32_t) i;

         for (j = qs_inf->small_primes; j < top && extra_bits < sieve[i]; j = jend)
         {
            jend = FLINT_MIN(j + QS_TRIAL_CHUNK, top);

            for (t = j; t < jend; t++)
            {
               uint32_t x1 = (i32 + p32[t] - (uint32_t) soln1[t]) * inv32[t];
               uint32_t x2 = (i32 + p32[t] - (uint32_t) soln2[t]) * inv32[t];
               hit[t - j] = (x1 <= lim32[t]) | (x2 <= lim32[t]);
            }

            for (t = j; t < jend; t++)
            {
               if (hit[t - j] && soln2[t] != 0)
               {
                  if (num_factors >= qs_inf->max_factors)
                     goto cleanup;

                  prime = factor_base[t].p;
                  exp = QS_REMOVE(prime, factor_base[t].pinv);
                  extra_bits += factor_base[t].size;
                  factor[num_factors].ind = t;
                  factor[num_factors++].exp = exp;
               }
            }
         }
      }

      /* the factors must be in order of their index */
      qsieve_sort_factors(factor, num_factors);

#if QS_DEBUG & 128
      if (num_factors > 0)
         flint_printf("\n");
#endif
      if (use_ll ? (rhi == 0 && rlo == 1)
                 : (fmpz_cmp_ui(res, 1) == 0 || fmpz_cmp_si(res, -1) == 0)) /* We've found a relation */
      {
#if QS_DEBUG
         if (qs_inf->full_relation % 100 == 0)
            flint_printf("%ld relations\n", qs_inf->full_relation);
#endif
         /* set sign amongst small factors */
         if (use_ll)
            small[2] = neg;
         else if (fmpz_cmp_si(res, -1) == 0)
            small[2] = 1;
         else
            small[2] = 0;

         poly->num_factors = num_factors;

#if FLINT_USES_PTHREAD
         pthread_mutex_lock(&qs_inf->mutex);
#endif
         qsieve_add_relation(qs_inf, 1, 1, Y, poly);

         qs_inf->full_relation++;

#if FLINT_USES_PTHREAD
         pthread_mutex_unlock(&qs_inf->mutex);
#endif
         relations++;
      } else /* not a relation, perhaps a partial? */
      {
          ulong r = 0;

          /* set sign */
          if (use_ll)
              small[2] = neg;
          else if (fmpz_sgn(res) < 0)
          {
              fmpz_neg(res, res);
              small[2] = 1;
          } else
              small[2] = 0;

          if (use_ll ? (rhi == 0) : fmpz_abs_fits_ui(res))
              r = use_ll ? rlo : fmpz_get_ui(res);

          /*
             The cofactor r has no prime factors in the factor base.  A large
             prime is taken heuristically to be below lp_bound, a small
             multiple of the largest factor base prime (so r is prime if it
             is below lp_bound).  Skip values not coprime with the
             multiplier, as these will lead to factors of kn, not n.
          */
          if (r != 0 && r < qs_inf->lp_bound && n_gcd(r, qs_inf->k) == 1)
          {
              poly->num_factors = num_factors;

#if FLINT_USES_PTHREAD
              pthread_mutex_lock(&qs_inf->mutex);
#endif
              /* store this partial */
              qsieve_add_relation(qs_inf, r, 1, Y, poly);
              qsieve_add_partial_edge(qs_inf, r, 1);

#if FLINT_USES_PTHREAD
              pthread_mutex_unlock(&qs_inf->mutex);
#endif
          }
      }

      goto cleanup;
   }

cleanup:
   fmpz_clear(X);
   fmpz_clear(Y);
   fmpz_clear(res);
   fmpz_clear(p);

   return relations;
}

#if FLINT64
#define QS_SCAN_MASK UWORD(0x8080808080808080)
#else
#define QS_SCAN_MASK UWORD(0x80808080)
#endif

/*
    scan sieve array for possible candidates for smooth relations and process
    those relations to check for smoothness
*/
slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
{
    slong i = 0, j = 0;
    ulong * sieve2 = (ulong *) sieve;
    unsigned char bits = qs_inf->sieve_bits;
    slong rels = 0;

    slong nwords = qs_inf->sieve_size / sizeof(ulong);

    while (j < nwords)
    {
        /*
           scan 4 or 8 bytes at once for sieve entries over the threshold,
           which have their top bit set (see qsieve_init_with_tune), and skip
           four words at a time while possible
        */
        while (j + 4 <= nwords && ((sieve2[j] | sieve2[j + 1] | sieve2[j + 2]
                                  | sieve2[j + 3]) & QS_SCAN_MASK) == 0)
        {
            j += 4;
        }

        while ((sieve2[j] & QS_SCAN_MASK) == 0)
        {
            j++; /* advance to next word */
        }

        i = j * sizeof(ulong);

        /* check bytes individually in word */
        while (i < (j + 1) * sizeof(ulong) && i < qs_inf->sieve_size)
        {
            /* if we are over the threshold, check candidate for smoothness */
            if (sieve[i] > bits)
               rels += qsieve_evaluate_candidate(qs_inf, i, sieve, poly);

            i++;
        }

        j++;
    }

    return rels;
}

/* procedure to call polynomial initialization and sieving procedure */

typedef struct
{
    qs_s * inf;
    unsigned char * sieve;
    slong thread_idx;
    qs_poly_s * thread_poly;
    unsigned char * thread_sieve;
    slong rels;
}
_worker_arg_struct;

static void qsieve_collect_relations_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    qs_s * qs_inf = arg->inf;
    qs_poly_s * thread_poly = arg->thread_poly;
    unsigned char * thread_sieve = arg->thread_sieve;
    slong j, iterations = (1 << (qs_inf->s - 1));
    /*
       Without other threads, sieve directly with the roots held in qs_inf,
       which are updated in place for each polynomial, instead of copying
       them (a copy of two words per factor base prime per polynomial).
    */
    int alias = (qs_inf->num_handles == 0);
    int * own_soln1 = thread_poly->soln1;
    int * own_soln2 = thread_poly->soln2;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&qs_inf->mutex);
#endif
        j = qs_inf->index_j;
        qs_inf->index_j = j + 1;
        if (j < iterations)
        {
           /* copy poly data for thread we are in */
            if (j > 0)
                qsieve_init_poly_next(qs_inf, j);

            if (alias)
            {
                fmpz_set(thread_poly->B, qs_inf->B);
                thread_poly->soln1 = qs_inf->soln1;
                thread_poly->soln2 = qs_inf->soln2;
            }
            else
                qsieve_poly_copy(thread_poly, qs_inf);
            thread_poly->C_valid = 0;
        }
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&qs_inf->mutex);
#endif

        if (j >= iterations)
        {
            thread_poly->soln1 = own_soln1;
            thread_poly->soln2 = own_soln2;
            return;
        }

        if (qs_inf->use_buckets)
           qsieve_do_sieving3(qs_inf, thread_sieve, thread_poly);
        else if (qs_inf->sieve_size < 2*QS_L1_SIZE)
           qsieve_do_sieving(qs_inf, thread_sieve, thread_poly);
        else
           qsieve_do_sieving2(qs_inf, thread_sieve, thread_poly);

        arg->rels += qsieve_evaluate_sieve(qs_inf, thread_sieve, thread_poly);
    }
}


slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)
{
    slong i;
    _worker_arg_struct * args;
    slong relations;
    const thread_pool_handle* handles = qs_inf->handles;
    slong num_handles = qs_inf->num_handles;

    args = (_worker_arg_struct *) flint_malloc((1 + num_handles)
                                                  *sizeof(_worker_arg_struct));

    qs_inf->index_j = 0;
    qsieve_init_poly_first(qs_inf);

    for (i = 0; i <= num_handles; i++)
    {
        args[i].inf = qs_inf;
        args[i].thread_idx = i;
        args[i].thread_poly = qs_inf->poly + i;
        args[i].thread_sieve = sieve + (qs_inf->sieve_size + sizeof(ulong) + 64)*i;
        args[i].rels = 0;
    }

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], 0,
                                    qsieve_collect_relations_worker, &args[i]);
    }

    qsieve_collect_relations_worker(&args[num_handles]);

    relations = args[num_handles].rels;

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
        relations += args[i].rels;
    }

    flint_free(args);

    return relations;
}
