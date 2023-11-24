/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "fmpz.h"
#include "qsieve.h"

#define HASH_MULT (2654435761U)       /* hash function, taken from 'msieve' */
#define HASH(a) ((ulong)((((unsigned int) a) * HASH_MULT) >> (12)))

/******************************************************************************
 *
 *  Some helper function, used for debugging
 *
 *****************************************************************************/

/*
    Display a relation for debugging purposes
*/
void qsieve_display_relation(qs_t qs_inf, relation_t a)
{
    slong i;

    flint_printf("%wu ", a.lp);

    for (i = 0; i < qs_inf->small_primes; i++)
        flint_printf("%wd ", a.small[i]);

    flint_printf("%wd ", a.num_factors);

    for (i = 0; i < a.num_factors; i++)
        flint_printf("%wd %wu ", a.factor[i].ind, a.factor[i].exp);

    fmpz_print(a.Y);
    flint_printf("\n");
}

/*
    Check a relation is valid (debugging)
*/
int qsieve_is_relation(qs_t qs_inf, relation_t a)
{
    slong i;
    fmpz_t temp, temp2;
    fmpz_init(temp);
    fmpz_init_set_ui(temp2, 1);

    for (i = 0; i < qs_inf->small_primes; i++)
    {
        fmpz_set_si(temp, qs_inf->factor_base[i].p);
        fmpz_pow_ui(temp, temp, a.small[i]);
        fmpz_mul(temp2, temp2, temp);
    }

    if (a.num_factors > qs_inf->max_factors)
    {
        return 0;
    }

    for (i = 0; i < a.num_factors; i++)
    {
        fmpz_set_ui(temp, qs_inf->factor_base[a.factor[i].ind].p);
        fmpz_pow_ui(temp, temp, a.factor[i].exp);
        fmpz_mul(temp2, temp2, temp);
    }

    fmpz_mul_ui(temp2, temp2, a.lp);
    fmpz_pow_ui(temp, a.Y, UWORD(2));
    fmpz_mod(temp, temp, qs_inf->kn);
    fmpz_mod(temp2, temp2, qs_inf->kn);

    if (fmpz_cmp(temp, temp2) != 0)
    {
        return 0;
    }

    fmpz_clear(temp);
    fmpz_clear(temp2);

    return 1;
}

/*
    Write partial or full relation to file

    The layout is as follows:
    total write size of relation (including this write size)
    large prime             (1 * mp_limb_t)
    number of small primes  (1 * slong)
    small primes            (number of small primes * slong)
    number of factors       (1 * slong)
    (factor, exponent)      (number of factors * fac_t)
    Y->_mp_size             (1 * slong)
    Y->_mp_d                (Y->_mp_size * mp_limb_t)
 */
void qsieve_write_to_file(qs_t qs_inf, mp_limb_t prime, const fmpz_t Y, const qs_poly_t poly)
{
    slong num_factors = poly->num_factors;
    slong * small = poly->small;
    fac_t * factor = poly->factor;
    slong Ysz;
    slong write_size;

    /* Get size of Y */
    Ysz = COEFF_IS_MPZ(*Y) ? COEFF_TO_PTR(*Y)->_mp_size : FLINT_SGN(*Y);

    /* Write size of relation */
    write_size =
        sizeof(slong)                           /* total write size */
        + sizeof(mp_limb_t)                     /* large prime */
        + sizeof(slong)                         /* number of small primes */
        + sizeof(slong) * qs_inf->small_primes  /* small primes */
        + sizeof(slong)                         /* number of factors */
        + sizeof(fac_t) * num_factors           /* factors */
        + sizeof(slong)                         /* Y->_mp_size */
        + sizeof(mp_limb_t) * (Ysz != 0 ? FLINT_ABS(Ysz) : 1); /* Y->_mp_d */
    QS_SIQS_FWRITE(&write_size, sizeof(slong), 1, qs_inf);

    /* Write large prime */
    QS_SIQS_FWRITE(&prime, sizeof(mp_limb_t), 1, qs_inf);

    /* NOTE: We do not have to write small primes. */
    /* Write number of small primes */
    QS_SIQS_FWRITE(&qs_inf->small_primes, sizeof(slong), 1, qs_inf);

    /* Write small primes */
    QS_SIQS_FWRITE(small, sizeof(slong), qs_inf->small_primes, qs_inf);

    /* Write number of factors */
    QS_SIQS_FWRITE(&num_factors, sizeof(slong), 1, qs_inf);

    /* Write factors and exponents */
    QS_SIQS_FWRITE(factor, sizeof(fac_t), num_factors, qs_inf);

    /* Write Y->_mp_size (or mock it) */
    QS_SIQS_FWRITE(&Ysz, sizeof(slong), 1, qs_inf);

    /* Write Y->_mp_d (or mock it) */
    if (!COEFF_IS_MPZ(*Y))
    {
        slong abslimb = FLINT_ABS(*Y);

        /* Write mock Y->_mp_d */
        QS_SIQS_FWRITE(&abslimb, sizeof(mp_limb_t), 1, qs_inf);
    }
    else
    {
        mp_srcptr Yd = COEFF_TO_PTR(*Y)->_mp_d;

        /* Write Y->_mp_d */
        QS_SIQS_FWRITE(Yd, sizeof(mp_limb_t), FLINT_ABS(Ysz), qs_inf);
    }
}

/******************************************************************************
 *
 *  Hash table
 *
 *****************************************************************************/

/*
   Hash table used to keep count of large primes, idea is taken from msieve
   Each new prime is filled at last unoccupied position in array and primes
   which have same hash value are linked with each other keeping offset
*/

/*
   return a pointer to location of 'prime' in table if it exists else
   create an entry for it and return pointer to that
*/
hash_t * qsieve_get_table_entry(qs_t qs_inf, mp_limb_t prime)
{
    mp_limb_t offset, first_offset;
    hash_t * entry;
    mp_limb_t * hash_table =  qs_inf->hash_table;
    hash_t * table = qs_inf->table;
    slong table_size = qs_inf->table_size;

    /* reallocate table if not large enough */
    if (3*qs_inf->vertices/2 + 1 >= table_size)
    {
        table_size *= 1.4;
        table = flint_realloc(table, table_size*sizeof(hash_t));
        qs_inf->table_size = table_size;
        qs_inf->table = table;
    }

    /* find first offset with that hash */
    first_offset = HASH(prime);
    offset = hash_table[first_offset];

    /* check linked offsets to see if prime is there, return if so */
    while (offset != 0)
    {
        entry = table + offset;
        if (entry->prime == prime)
            break;
        offset = entry->next;
    }

    /* if we didn't find it, make a new entry in hash table and return it */
    if (offset == 0)
    {
        qs_inf->vertices++;
        entry = table + qs_inf->vertices;
        entry->prime = prime;
        entry->next = hash_table[first_offset];
        entry->count = 0;
        hash_table[first_offset] = qs_inf->vertices;
    }

    return entry;
}

/*
   add prime to hashtable, increase size of table if necessary
   and increment count for the added prime
*/
void qsieve_add_to_hashtable(qs_t qs_inf, mp_limb_t prime)
{
    hash_t * entry;

    entry = qsieve_get_table_entry(qs_inf, prime);
    entry->count++;
}

/******************************************************************************
 *
 *  Large prime functionality
 *
 *****************************************************************************/

/*
   given a string representing a relation, parse it to
   obtain relation
*/
relation_t qsieve_parse_relation(qs_t qs_inf)
{
    relation_t rel;
    slong Ysz;

    /* NOTE: write_size and large prime is already read in
     * qsieve_process_relation. */

    /* Get large prime (is always one) */
    rel.lp = UWORD(1);

    /* NOTE: We can use qs_inf->small_primes here instead of reading. */
    /* Get number of small primes */
    QS_SIQS_FREAD_NORES(&rel.small_primes, sizeof(slong), 1, qs_inf);

    /* Get small primes */
    rel.small = flint_malloc(rel.small_primes * sizeof(slong));
    QS_SIQS_FREAD_NORES(rel.small, sizeof(slong), rel.small_primes, qs_inf);

    /* Get number of factors */
    QS_SIQS_FREAD_NORES(&rel.num_factors, sizeof(slong), 1, qs_inf);

    /* Get factors */
    rel.factor = flint_malloc(rel.num_factors * sizeof(fac_t));
    QS_SIQS_FREAD_NORES(rel.factor, sizeof(fac_t), rel.num_factors, qs_inf);

    /* Get Ysz */
    Ysz = 0;
    QS_SIQS_FREAD_NORES(&Ysz, sizeof(slong), 1, qs_inf);

    /* Get Y */
    fmpz_init(rel.Y);
    if (FLINT_ABS(Ysz) <= 1)
    {
        mp_limb_t abslimb = 0;

        QS_SIQS_FREAD_NORES(&abslimb, sizeof(mp_limb_t), 1, qs_inf);

#if COEFF_MAX != -COEFF_MIN
# error
#endif
        fmpz_set_ui(rel.Y, abslimb);
        if (Ysz < 0)
            fmpz_neg(rel.Y, rel.Y);
    }
    else
    {
        mpz_ptr mY = _fmpz_new_mpz();

        mY->_mp_size = Ysz;

        if (mY->_mp_alloc < FLINT_ABS(Ysz))
            _mpz_realloc(mY, FLINT_ABS(Ysz));

        QS_SIQS_FREAD_NORES(mY->_mp_d, sizeof(mp_limb_t), FLINT_ABS(Ysz), qs_inf);
        *rel.Y = PTR_TO_COEFF(mY);
    }

    return rel;
}

/*
   given two partials with same large prime, merge them to
   obtain a full relation
*/
relation_t qsieve_merge_relation(qs_t qs_inf, relation_t a, relation_t b)
{
    slong i = 0, j = 0, k = 0;
    relation_t  c;
    fmpz_t temp;

    c.lp = UWORD(1);
    c.small = flint_malloc(qs_inf->small_primes * sizeof(slong));
    c.factor = flint_malloc(qs_inf->max_factors * sizeof(fac_t));
    fmpz_init(c.Y);

    for (i = 0; i < qs_inf->small_primes; i++)
        c.small[i] = (a.small[i] + b.small[i]);

    i = 0;

    while (i < a.num_factors && j < b.num_factors)
    {
        if (a.factor[i].ind == b.factor[j].ind)
        {
            c.factor[k].ind = a.factor[i].ind;
            c.factor[k++].exp = a.factor[i++].exp + b.factor[j++].exp;
        }
        else if (a.factor[i].ind < b.factor[j].ind)
        {
            c.factor[k].ind = a.factor[i].ind;
            c.factor[k++].exp = a.factor[i++].exp;
        }
        else
        {
           c.factor[k].ind = b.factor[j].ind;
           c.factor[k++].exp = b.factor[j++].exp;
        }

        if (k >= qs_inf->max_factors)
        {
            flint_printf("more than max_factor !!\n");
            flint_abort();
        }
    }

    while (i < a.num_factors)
    {
        c.factor[k].ind = a.factor[i].ind;
        c.factor[k++].exp = a.factor[i++].exp;

        if (k >= qs_inf->max_factors)
        {
            flint_printf("more than max_factor !!\n");
            flint_abort();
        }
    }

    while (j < b.num_factors)
    {
        c.factor[k].ind = b.factor[j].ind;
        c.factor[k++].exp = b.factor[j++].exp;

        if (k >= qs_inf->max_factors)
        {
            flint_printf("more than max_factor !!\n");
            flint_abort();
        }
    }

    c.num_factors = k;
    c.small_primes = qs_inf->small_primes;

    fmpz_init_set_ui(temp, a.lp);

    if (fmpz_invmod(temp, temp, qs_inf->kn) == 0)
    {
        flint_printf("Inverse doesn't exist !!\n");
        flint_abort();
    }

    fmpz_mul(c.Y, a.Y, b.Y);
    fmpz_mul(c.Y, c.Y, temp);
    if (fmpz_cmp(qs_inf->kn, c.Y) <= 0)
        fmpz_mod(c.Y, c.Y, qs_inf->kn);
    fmpz_clear(temp);

    return c;
}

/*
   compare two relations in the following order,
   large_prime, number of factors, factor, small_prime
*/
int qsieve_compare_relation(const void * a, const void * b)
{
    slong i;
    relation_t * r1 = (relation_t *) a;
    relation_t * r2 = (relation_t *) b;

    if (r1->lp > r2->lp)
        return 1;

    if (r1->lp < r2->lp)
        return -1;

    if (r1->num_factors > r2->num_factors)
        return 1;

    if (r1->num_factors < r2->num_factors)
        return -1;

    for (i = 0; i < r1->num_factors; i++)
    {
        if (r1->factor[i].ind > r2->factor[i].ind)
            return 1;

        if (r1->factor[i].ind < r2->factor[i].ind)
            return -1;

        if (r1->factor[i].exp > r2->factor[i].exp)
            return 1;

        if (r1->factor[i].exp < r2->factor[i].exp)
            return -1;
    }

    for (i = 0; i < r1->small_primes; i++)
    {
        if (r1->small[i] > r2->small[i])
            return 1;

        if (r1->small[i] < r2->small[i])
            return -1;
    }

    return 0;
}

/*
   given a list of relations, remove duplicate relations from it
*/
int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations)
{
    slong i, j;

    if (num_relations < 2)
        return 1;

    qsort(rel_list, (size_t) num_relations, sizeof(relation_t), qsieve_compare_relation);

    for (i = 1, j = 0; i < num_relations; i++)
    {
        if (qsieve_compare_relation(rel_list + j, rel_list + i) == 0)
        {
            rel_list[i].num_factors = 0;
            flint_free(rel_list[i].small);
            flint_free(rel_list[i].factor);
            fmpz_clear(rel_list[i].Y);
        } else
        {
            rel_list[++j] = rel_list[i];
        }
    }

    j++;

#if QS_DEBUG
    flint_printf("%wd duplicates out of %wd\n", num_relations - j, num_relations);
#endif

    return j;
}

/*
   give a list of relations, add those relations to matrix
*/
void qsieve_insert_relation(qs_t qs_inf, relation_t * rel_list, slong num_relations)
{
    slong i, j, num_factors, fac_num;
    slong * small;
    slong * curr_rel;
    fac_t * factor;
    la_col_t * matrix = qs_inf->matrix;

    qs_inf->num_relations = 0;

    for (j = 0; j < num_relations; j++)
    {
        small = rel_list[j].small;
        num_factors = rel_list[j].num_factors;
        factor = rel_list[j].factor;
        curr_rel = qs_inf->curr_rel;
        fac_num = 0;

        clear_col(matrix + j);

        for (i = 0; i < qs_inf->small_primes; i++)
        {
            if (small[i] & 1) insert_col_entry(matrix + j, i);

            if (small[i])
            {
                curr_rel[2*fac_num + 1] = i;
                curr_rel[2*fac_num + 2] = small[i];
                fac_num++;
            }
        }

        for (i = 0; i < num_factors; i++)
        {
            if (factor[i].exp & 1) insert_col_entry(matrix + j, factor[i].ind);
            curr_rel[2*fac_num + 1] = factor[i].ind;
            curr_rel[2*fac_num + 2] = factor[i].exp;
            fac_num++;
        }

        curr_rel[0] = fac_num;

        matrix[j].orig = qs_inf->num_relations;

        fmpz_set(qs_inf->Y_arr + qs_inf->num_relations, rel_list[j].Y);

        qs_inf->curr_rel += qs_inf->max_factors*2;
        qs_inf->num_relations++;
    }

    qs_inf->columns = qs_inf->num_relations;
}

/*
   process relations from the file
*/
int qsieve_process_relation(qs_t qs_inf)
{
    slong i, num_relations = 0, num_relations2;
    slong rel_list_length;
    slong rlist_length;
    mp_limb_t prime;
    hash_t * entry;
    mp_limb_t * hash_table = qs_inf->hash_table;
    slong rel_size = 50000;
    relation_t * rel_list = (relation_t *) flint_malloc(rel_size * sizeof(relation_t));
    relation_t * rlist;
    int done = 0;

    QS_SIQS_FCLOSE(qs_inf);
    QS_SIQS_FOPEN_R(qs_inf);

#if QS_DEBUG & 64
    flint_printf("Getting relations\n");
#endif

    while (1)
    {
        int read_size;
        slong relation_size = 0;

        QS_SIQS_FREAD(read_size, &relation_size, sizeof(slong), 1, qs_inf);

        if (read_size != sizeof(slong))
            break;

        QS_SIQS_FREAD_NORES(&prime, sizeof(mp_limb_t), 1, qs_inf);
        entry = qsieve_get_table_entry(qs_inf, prime);

        if (num_relations == rel_size)
        {
           rel_list = (relation_t *) flint_realloc(rel_list, 2 * rel_size * sizeof(relation_t));
           rel_size *= 2;
        }

        if (prime == 1 || entry->count >= 2)
        {
            rel_list[num_relations] = qsieve_parse_relation(qs_inf);
            rel_list[num_relations].lp = prime;
            num_relations++;
        }
        else
        {
            /* We have to get to the next relation in the file. We have already
             * read write_size (is a slong) and large prime (is an mp_limb_t).*/
            QS_SIQS_FSEEK_SEEK_CUR(qs_inf, relation_size - sizeof(slong) - sizeof(mp_limb_t));
        }
    }

    QS_SIQS_FCLOSE(qs_inf);

#if QS_DEBUG & 64
    flint_printf("Removing duplicates\n");
#endif

    num_relations = qsieve_remove_duplicates(rel_list, num_relations);
    rel_list_length = num_relations;

#if QS_DEBUG & 64
    flint_printf("Merging relations\n");
#endif

    rlist = flint_malloc(num_relations * sizeof(relation_t));
    memset(hash_table, 0, (1 << 20) * sizeof(mp_limb_t));
    qs_inf->vertices = 0;

    rlist_length = 0;
    for (i = 0; i < num_relations; i++)
    {
        if (rel_list[i].lp == UWORD(1))
        {
            rlist[rlist_length++] = rel_list[i];
        }
        else
        {
            entry = qsieve_get_table_entry(qs_inf, rel_list[i].lp);

            if (entry->count == 0) entry->count = i;
            else
            {
                if (fmpz_fdiv_ui(qs_inf->kn, rel_list[i].lp) == 0)
                {
                   qs_inf->small_factor = rel_list[i].lp;

                   done = -1;
                   goto cleanup;
                }
                rlist[rlist_length++] = qsieve_merge_relation(qs_inf, rel_list[i], rel_list[entry->count]);
            }
        }
    }

    num_relations = rlist_length;

#if QS_DEBUG & 64
    flint_printf("Sorting relations\n");
#endif

    if (rlist_length < qs_inf->num_primes + qs_inf->ks_primes + qs_inf->extra_rels)
    {
       qs_inf->edges -= 100;
       done = 0;
       QS_SIQS_FCLOSE(qs_inf);
       QS_SIQS_FOPEN_A(qs_inf);
    } else
    {
       done = 1;
       num_relations2 = qs_inf->num_primes + qs_inf->ks_primes + qs_inf->extra_rels;
       qsort(rlist, (size_t) num_relations2, sizeof(relation_t), qsieve_compare_relation);
       qsieve_insert_relation(qs_inf, rlist, num_relations2);
    }

cleanup:

    for (i = 0; i < rel_list_length; i++)
    {
        /* it looks like rlist stole our data if rel_list[i].lp == UWORD(1)) */
        if (rel_list[i].lp != UWORD(1))
        {
            flint_free(rel_list[i].small);
            flint_free(rel_list[i].factor);
            fmpz_clear(rel_list[i].Y);
        }
    }
    flint_free(rel_list);

    for (i = 0; i < rlist_length; i++)
    {
       flint_free(rlist[i].small);
       flint_free(rlist[i].factor);
       fmpz_clear(rlist[i].Y);
    }
    flint_free(rlist);

    return done;
}

