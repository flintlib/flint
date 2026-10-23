/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "qsieve.h"

#define HASH_MULT (2654435761U)       /* hash function, taken from 'msieve' */
#define HASH(a) ((ulong)((((unsigned int) a) * HASH_MULT) >> qs_inf->hash_shift))

/******************************************************************************
 *
 *  Relation store
 *
 *****************************************************************************/

/*
   The relations found by sieving are kept in memory until the linear algebra
   consumes them.  See the comment on the store in qsieve.h for the layout.
*/

#define QS_REL_ALLOC_MIN 256          /* initial number of relations */
#define QS_REL_FACTOR_ALLOC_MIN 4096  /* initial number of factor entries */

void qsieve_relations_init(qs_t qs_inf)
{
    qs_inf->rel_num = 0;
    qs_inf->rel_alloc = QS_REL_ALLOC_MIN;
    qs_inf->relations = flint_malloc(qs_inf->rel_alloc*sizeof(qs_rel_s));
    qs_inf->rel_small = flint_malloc(qs_inf->rel_alloc*qs_inf->small_primes);

    qs_inf->rel_factor_num = 0;
    qs_inf->rel_factor_alloc = QS_REL_FACTOR_ALLOC_MIN;
    qs_inf->rel_factor = flint_malloc(qs_inf->rel_factor_alloc*sizeof(uint32_t));
}

/* discard all relations found so far, keeping the allocations */
void qsieve_relations_reset(qs_t qs_inf)
{
    slong i;

    for (i = 0; i < qs_inf->rel_num; i++)
        fmpz_clear(&qs_inf->relations[i].Y);

    qs_inf->rel_num = 0;
    qs_inf->rel_factor_num = 0;
}

void qsieve_relations_clear(qs_t qs_inf)
{
    qsieve_relations_reset(qs_inf);

    flint_free(qs_inf->relations);
    flint_free(qs_inf->rel_small);
    flint_free(qs_inf->rel_factor);

    qs_inf->relations = NULL;
    qs_inf->rel_small = NULL;
    qs_inf->rel_factor = NULL;
    qs_inf->rel_alloc = 0;
    qs_inf->rel_factor_alloc = 0;
}

/*
   Add a full (prime == prime2 == 1) or partial relation to the store.  A
   partial with a single large prime has prime2 == 1.  (The relations, and
   the processing in qsieve_process_relation, allow for two large primes,
   but the sieve currently only produces partials with a single one.)

   Called from the sieving threads, so the caller must hold qs_inf->mutex.
*/
void qsieve_add_relation(qs_t qs_inf, ulong prime, ulong prime2, const fmpz_t Y, const qs_poly_t poly)
{
    slong num_factors = poly->num_factors;
    slong small_primes = qs_inf->small_primes;
    qs_rel_s * rel;

    if (qs_inf->rel_num == qs_inf->rel_alloc)
    {
        slong alloc = 2*qs_inf->rel_alloc;

        qs_inf->relations = flint_realloc(qs_inf->relations,
                                          alloc*sizeof(qs_rel_s));
        qs_inf->rel_small = flint_realloc(qs_inf->rel_small,
                                          alloc*small_primes);
        qs_inf->rel_alloc = alloc;
    }

    if (qs_inf->rel_factor_num + num_factors > qs_inf->rel_factor_alloc)
    {
        slong alloc = FLINT_MAX(2*qs_inf->rel_factor_alloc,
                                qs_inf->rel_factor_num + num_factors);

        qs_inf->rel_factor = flint_realloc(qs_inf->rel_factor,
                                           alloc*sizeof(uint32_t));
        qs_inf->rel_factor_alloc = alloc;
    }

    rel = qs_inf->relations + qs_inf->rel_num;

    rel->lp = prime;
    rel->lp2 = prime2;
    rel->num_factors = num_factors;
    rel->factor_offset = qs_inf->rel_factor_num;
    fmpz_init_set(&rel->Y, Y);

    {
        slong i;
        unsigned char * s = qs_inf->rel_small + qs_inf->rel_num*small_primes;

        for (i = 0; i < small_primes; i++)
        {
            FLINT_ASSERT(poly->small[i] >= 0 && poly->small[i] < 256);
            s[i] = poly->small[i];
        }
    }
    {
        slong i;
        uint32_t * f = qs_inf->rel_factor + qs_inf->rel_factor_num;

        for (i = 0; i < num_factors; i++)
        {
            FLINT_ASSERT(poly->factor[i].ind < (WORD(1) << 24));
            FLINT_ASSERT(poly->factor[i].exp > 0 && poly->factor[i].exp < 256);
            f[i] = (((uint32_t) poly->factor[i].ind) << 8) | poly->factor[i].exp;
        }
    }

    qs_inf->rel_factor_num += num_factors;
    qs_inf->rel_num++;
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
hash_t * qsieve_get_table_entry(qs_t qs_inf, ulong prime)
{
    ulong offset, first_offset;
    hash_t * entry;
    ulong * hash_table =  qs_inf->hash_table;
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
        entry->parent = qs_inf->vertices;
        hash_table[first_offset] = qs_inf->vertices;
    }

    return entry;
}

/*
   add prime to hashtable, increase size of table if necessary
   and increment count for the added prime
*/
void qsieve_add_to_hashtable(qs_t qs_inf, ulong prime)
{
    hash_t * entry;

    entry = qsieve_get_table_entry(qs_inf, prime);
    entry->count++;
}

/* root of x in the union-find forest, with path halving */
static ulong
qsieve_uf_find(hash_t * table, ulong x)
{
    while (table[x].parent != x)
    {
        table[x].parent = table[table[x].parent].parent;
        x = table[x].parent;
    }

    return x;
}

/*
   Record the edge prime -- prime2 of the large prime graph for a partial
   relation (prime2 = 1 for a single large prime).  The partial completes a
   new independent cycle, and so contributes a full relation, exactly when
   its two endpoints are already connected; this is tracked with a union-find
   structure on the vertices.

   Called from the sieving threads, so the caller must hold qs_inf->mutex.
*/
void qsieve_add_partial_edge(qs_t qs_inf, ulong prime, ulong prime2)
{
    ulong u, v;
    hash_t * table;

    /* the second lookup may reallocate the table, so keep indices */
    u = qsieve_get_table_entry(qs_inf, prime) - qs_inf->table;
    v = qsieve_get_table_entry(qs_inf, prime2) - qs_inf->table;
    table = qs_inf->table;

    table[u].count++;
    table[v].count++;
    qs_inf->edges++;

    u = qsieve_uf_find(table, u);
    v = qsieve_uf_find(table, v);

    if (u == v)
        qs_inf->num_cycles++;
    else
        table[u].parent = v;
}

/******************************************************************************
 *
 *  Large prime functionality
 *
 *****************************************************************************/

/*
   take a copy of the i-th relation in the store

   The store keeps ownership of its own copy: the relation returned here is
   handed to the merging and filtering code below, which consumes and frees it.
*/
relation_t qsieve_get_relation(qs_t qs_inf, slong i)
{
    qs_rel_s * stored = qs_inf->relations + i;
    relation_t rel;

    rel.lp = stored->lp;
    rel.lp2 = stored->lp2;
    rel.small_primes = qs_inf->small_primes;
    rel.num_factors = stored->num_factors;

    rel.small = flint_malloc(rel.small_primes * sizeof(slong));
    {
        slong j;
        const unsigned char * s = qs_inf->rel_small + i*qs_inf->small_primes;

        for (j = 0; j < rel.small_primes; j++)
            rel.small[j] = s[j];
    }

    rel.factor = flint_malloc(rel.num_factors * sizeof(fac_t));
    {
        slong j;
        const uint32_t * f = qs_inf->rel_factor + stored->factor_offset;

        for (j = 0; j < rel.num_factors; j++)
        {
            rel.factor[j].ind = f[j] >> 8;
            rel.factor[j].exp = f[j] & 255;
        }
    }

    fmpz_init_set(rel.Y, &stored->Y);

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
    c.lp2 = UWORD(1);
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
            flint_throw(FLINT_ERROR, "more than max_factor !!\n");
        }
    }

    while (i < a.num_factors)
    {
        c.factor[k].ind = a.factor[i].ind;
        c.factor[k++].exp = a.factor[i++].exp;

        if (k >= qs_inf->max_factors)
        {
            flint_throw(FLINT_ERROR, "more than max_factor !!\n");
        }
    }

    while (j < b.num_factors)
    {
        c.factor[k].ind = b.factor[j].ind;
        c.factor[k++].exp = b.factor[j++].exp;

        if (k >= qs_inf->max_factors)
        {
            flint_throw(FLINT_ERROR, "more than max_factor !!\n");
        }
    }

    c.num_factors = k;
    c.small_primes = qs_inf->small_primes;

    fmpz_init_set_ui(temp, a.lp);

    if (fmpz_invmod(temp, temp, qs_inf->kn) == 0)
    {
        flint_throw(FLINT_ERROR, "Inverse doesn't exist !!\n");
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

    if (r1->lp2 > r2->lp2)
        return 1;

    if (r1->lp2 < r2->lp2)
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

   The factorisations are stored one after the other in qs_inf->relation,
   each as a count followed by (index, exponent) pairs, with relation i
   starting at qs_inf->rel_offset[i].  (They used to be stored at a fixed
   stride of 2*max_factors, which silently overflowed for relations with
   more factors, and wasted memory for the others.)
*/
void qsieve_insert_relation(qs_t qs_inf, relation_t * rel_list, slong num_relations)
{
    slong i, j, num_factors, fac_num, need;
    slong * small;
    slong * curr_rel;
    fac_t * factor;
    la_col_t * matrix = qs_inf->matrix;

    qs_inf->num_relations = 0;

    need = 0;
    for (j = 0; j < num_relations; j++)
        need += 1 + 2*(qs_inf->small_primes + rel_list[j].num_factors);

    if (need > qs_inf->relation_alloc)
    {
        flint_free(qs_inf->relation);
        qs_inf->relation = flint_malloc(need*sizeof(slong));
        qs_inf->relation_alloc = need;
    }

    qs_inf->rel_offset = flint_realloc(qs_inf->rel_offset,
                                       FLINT_MAX(num_relations, 1)*sizeof(slong));

    curr_rel = qs_inf->relation;

    for (j = 0; j < num_relations; j++)
    {
        small = rel_list[j].small;
        num_factors = rel_list[j].num_factors;
        factor = rel_list[j].factor;
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

        qs_inf->rel_offset[qs_inf->num_relations] = curr_rel - qs_inf->relation;
        curr_rel += 2*fac_num + 1;
        qs_inf->num_relations++;
    }

    qs_inf->columns = qs_inf->num_relations;
}

static void
qsieve_relation_clear(relation_t * rel)
{
    flint_free(rel->small);
    flint_free(rel->factor);
    fmpz_clear(rel->Y);
}

static int
qsieve_fac_cmp(const void * a, const void * b)
{
    slong x = ((const fac_t *) a)->ind, y = ((const fac_t *) b)->ind;
    return (x > y) - (x < y);
}

static int
qsieve_ulong_cmp(const void * a, const void * b)
{
    ulong x = *(const ulong *) a, y = *(const ulong *) b;
    return (x > y) - (x < y);
}

/* index of p in the sorted array v of length n, which must contain it */
static slong
qsieve_vertex_index(const ulong * v, slong n, ulong p)
{
    slong lo = 0, hi = n - 1;

    while (lo < hi)
    {
        slong mid = lo + (hi - lo)/2;

        if (v[mid] < p)
            lo = mid + 1;
        else
            hi = mid;
    }

    return lo;
}

/*
   Combine the relations rels[L[0]], ..., rels[L[len - 1]] of a cycle of the
   large prime graph into a full relation.  The large primes of the cycle are
   W[0], ..., W[wlen - 1] (entries 1 are ignored), and each of them occurs in
   exactly two of the relations, so the product of the relations is a full
   relation times the square of prod W: its exponents are the sums of the
   exponents, and its square root is Y = prod Y_i / prod W mod kn.  Returns
   0, and sets qs_inf->small_factor, if one of the large primes is not
   coprime to kn.
*/
static int
qsieve_combine_relations(relation_t * res, qs_t qs_inf, const relation_t * rels,
                const slong * L, slong len, const ulong * W, slong wlen)
{
    slong i, j, k, total;
    fac_t * fac;
    fmpz_t t, w;

    total = 0;
    for (i = 0; i < len; i++)
        total += rels[L[i]].num_factors;

    fac = flint_malloc(FLINT_MAX(total, 1)*sizeof(fac_t));

    for (i = k = 0; i < len; i++)
    {
        memcpy(fac + k, rels[L[i]].factor, rels[L[i]].num_factors*sizeof(fac_t));
        k += rels[L[i]].num_factors;
    }

    qsort(fac, total, sizeof(fac_t), qsieve_fac_cmp);

    /* merge entries with equal index */
    for (i = 0, j = -1; i < total; i++)
    {
        if (j >= 0 && fac[j].ind == fac[i].ind)
            fac[j].exp += fac[i].exp;
        else
            fac[++j] = fac[i];
    }

    res->lp = UWORD(1);
    res->lp2 = UWORD(1);
    res->small_primes = qs_inf->small_primes;
    res->num_factors = j + 1;
    res->factor = fac;
    res->small = flint_calloc(qs_inf->small_primes, sizeof(slong));

    for (i = 0; i < len; i++)
        for (j = 0; j < qs_inf->small_primes; j++)
            res->small[j] += rels[L[i]].small[j];

    fmpz_init(res->Y);
    fmpz_init(t);
    fmpz_init_set_ui(w, 1);

    fmpz_one(res->Y);
    for (i = 0; i < len; i++)
    {
        fmpz_mul(res->Y, res->Y, rels[L[i]].Y);
        fmpz_mod(res->Y, res->Y, qs_inf->kn);
    }

    for (i = 0; i < wlen; i++)
    {
        if (W[i] != 1)
        {
            fmpz_mul_ui(w, w, W[i]);
            fmpz_mod(w, w, qs_inf->kn);
        }
    }

    if (!fmpz_invmod(t, w, qs_inf->kn))
    {
        /* some large prime divides kn */
        for (i = 0; i < wlen; i++)
        {
            if (W[i] != 1 && fmpz_fdiv_ui(qs_inf->kn, W[i]) == 0)
            {
                qs_inf->small_factor = W[i];
                break;
            }
        }

        fmpz_clear(t);
        fmpz_clear(w);
        qsieve_relation_clear(res);
        return 0;
    }

    fmpz_mul(res->Y, res->Y, t);
    fmpz_mod(res->Y, res->Y, qs_inf->kn);

    fmpz_clear(t);
    fmpz_clear(w);

    return 1;
}

/* cycles longer than this give denser columns and are not used */
#define QS_MAX_CYCLE_LENGTH 12

/*
   process the relations found so far

   The partial relations are the edges of a graph whose vertices are the
   large primes and 1: a single partial with large prime p is an edge 1 -- p,
   and a double partial with large primes p and q an edge p -- q.  Every
   cycle in this graph combines to a full relation.  We find a spanning
   forest by breadth first search, and every edge not in it closes a cycle
   with the tree paths from its endpoints, which gives a basis of the cycle
   space.  With single large primes only, this simply combines every partial
   with the first one having the same large prime.
*/
int qsieve_process_relation(qs_t qs_inf)
{
    slong i, j, num_relations, num_relations2, num_fulls, ne, nv, rlist_length;
    slong rlist_alloc, head, tail;
    relation_t * rel_list, * rlist;
    ulong * verts;
    slong * eu, * ev, * deg, * off, * adj, * queue, * par_v, * par_e, * depth;
    char * alive, * tree;
    slong L[2*QS_MAX_CYCLE_LENGTH + 2];
    ulong W[2*QS_MAX_CYCLE_LENGTH + 2];
    slong one = -1;
    int done = 0;

    /*
       Take the full relations, and the partials whose large primes all
       occur in at least two partials, since a vertex of degree one is not
       on any cycle.
    */
    rel_list = flint_malloc(FLINT_MAX(qs_inf->rel_num, 1)*sizeof(relation_t));
    num_relations = 0;

    for (i = 0; i < qs_inf->rel_num; i++)
    {
        ulong lp = qs_inf->relations[i].lp, lp2 = qs_inf->relations[i].lp2;

        if (lp != 1 && qsieve_get_table_entry(qs_inf, lp)->count < 2)
            continue;
        if (lp2 != 1 && qsieve_get_table_entry(qs_inf, lp2)->count < 2)
            continue;

        rel_list[num_relations++] = qsieve_get_relation(qs_inf, i);
    }

    /* sorts by large primes, so the full relations come first */
    num_relations = qsieve_remove_duplicates(rel_list, num_relations);

    for (num_fulls = 0; num_fulls < num_relations; num_fulls++)
        if (rel_list[num_fulls].lp != 1)
            break;

    /* the graph on the remaining (partial) relations */
    ne = num_relations - num_fulls;

    verts = flint_malloc(FLINT_MAX(2*ne + 1, 1)*sizeof(ulong));
    nv = 0;
    verts[nv++] = 1;
    for (i = 0; i < ne; i++)
    {
        verts[nv++] = rel_list[num_fulls + i].lp;
        verts[nv++] = rel_list[num_fulls + i].lp2;
    }

    qsort(verts, nv, sizeof(ulong), qsieve_ulong_cmp);

    for (i = j = 0; i < nv; i++)
        if (j == 0 || verts[j - 1] != verts[i])
            verts[j++] = verts[i];
    nv = j;
    one = 0; /* verts[0] == 1 */

    eu = flint_malloc(FLINT_MAX(ne, 1)*sizeof(slong));
    ev = flint_malloc(FLINT_MAX(ne, 1)*sizeof(slong));
    alive = flint_malloc(FLINT_MAX(ne, 1));
    tree = flint_calloc(FLINT_MAX(ne, 1), 1);
    deg = flint_calloc(nv, sizeof(slong));
    off = flint_calloc(nv + 1, sizeof(slong));
    adj = flint_malloc(FLINT_MAX(2*ne, 1)*sizeof(slong));
    queue = flint_malloc(nv*sizeof(slong));
    par_v = flint_malloc(nv*sizeof(slong));
    par_e = flint_malloc(nv*sizeof(slong));
    depth = flint_malloc(nv*sizeof(slong));

    for (i = 0; i < ne; i++)
    {
        eu[i] = qsieve_vertex_index(verts, nv, rel_list[num_fulls + i].lp);
        ev[i] = qsieve_vertex_index(verts, nv, rel_list[num_fulls + i].lp2);
        alive[i] = 1;
        deg[eu[i]]++;
        deg[ev[i]]++;
    }

    for (i = 0; i < ne; i++)
    {
        off[eu[i] + 1]++;
        off[ev[i] + 1]++;
    }

    for (i = 0; i < nv; i++)
        off[i + 1] += off[i];

    {
        slong * pos = flint_malloc(nv*sizeof(slong));

        for (i = 0; i < nv; i++)
            pos[i] = off[i];

        for (i = 0; i < ne; i++)
        {
            adj[pos[eu[i]]++] = i;
            adj[pos[ev[i]]++] = i;
        }

        flint_free(pos);
    }

    /* repeatedly remove vertices of degree one with their edge */
    head = tail = 0;
    for (i = 0; i < nv; i++)
        if (deg[i] == 1)
            queue[tail++] = i;

    while (head < tail)
    {
        slong x = queue[head++];

        if (deg[x] != 1)
            continue;

        for (j = off[x]; j < off[x + 1]; j++)
        {
            slong e = adj[j];

            if (alive[e])
            {
                slong y = (eu[e] == x) ? ev[e] : eu[e];

                alive[e] = 0;
                deg[x]--;
                deg[y]--;

                if (deg[y] == 1)
                    queue[tail++] = y;
                break;
            }
        }
    }

    /* breadth first search, starting from the vertex 1 */
    for (i = 0; i < nv; i++)
        par_v[i] = -2; /* unvisited */

    for (i = -1; i < nv; i++)
    {
        slong root = (i == -1) ? one : i;

        if (par_v[root] != -2 || deg[root] == 0)
            continue;

        head = tail = 0;
        queue[tail++] = root;
        par_v[root] = -1;
        par_e[root] = -1;
        depth[root] = 0;

        while (head < tail)
        {
            slong x = queue[head++];

            for (j = off[x]; j < off[x + 1]; j++)
            {
                slong e = adj[j];
                slong y = (eu[e] == x) ? ev[e] : eu[e];

                if (alive[e] && par_v[y] == -2)
                {
                    par_v[y] = x;
                    par_e[y] = e;
                    depth[y] = depth[x] + 1;
                    tree[e] = 1;
                    queue[tail++] = y;
                }
            }
        }
    }

    /* the full relations, then one relation for every non tree edge */
    rlist_alloc = num_fulls + ne + 1;
    rlist = flint_malloc(rlist_alloc*sizeof(relation_t));

    for (i = 0; i < num_fulls; i++)
        rlist[i] = rel_list[i];
    rlist_length = num_fulls;

    for (i = 0; i < ne && done == 0; i++)
    {
        slong x, y, len = 0, wlen = 0;

        if (!alive[i] || tree[i])
            continue;

        L[len++] = num_fulls + i;
        x = eu[i];
        y = ev[i];

        while (depth[x] > depth[y] && len <= QS_MAX_CYCLE_LENGTH)
        {
            W[wlen++] = verts[x];
            L[len++] = num_fulls + par_e[x];
            x = par_v[x];
        }

        while (depth[y] > depth[x] && len <= QS_MAX_CYCLE_LENGTH)
        {
            W[wlen++] = verts[y];
            L[len++] = num_fulls + par_e[y];
            y = par_v[y];
        }

        while (x != y && len <= QS_MAX_CYCLE_LENGTH)
        {
            W[wlen++] = verts[x];
            W[wlen++] = verts[y];
            L[len++] = num_fulls + par_e[x];
            L[len++] = num_fulls + par_e[y];
            x = par_v[x];
            y = par_v[y];
        }

        if (len > QS_MAX_CYCLE_LENGTH)
            continue;

        W[wlen++] = verts[x]; /* the common ancestor */

        if (!qsieve_combine_relations(rlist + rlist_length, qs_inf,
                                       rel_list, L, len, W, wlen))
        {
            done = -1;
            break;
        }

        rlist_length++;
    }


    if (done == 0)
    {
        /*
           Only the first num_primes factor base primes are sieved with (the
           factor base is computed with ks_primes more), so only those can
           appear in a relation, and extra_rels more relations than that are
           needed.
        */
        num_relations2 = qs_inf->num_primes + qs_inf->extra_rels;

        if (rlist_length < num_relations2)
        {
           /*
              Not enough (the estimate counts duplicates and overlong
              cycles): keep the relations we have and go back to sieving.
              Lower the estimate by the shortfall, so that we sieve for at
              least that many more relations before trying again.
           */
           qs_inf->num_cycles -= num_relations2 - rlist_length;
        }
        else
        {
           done = 1;
           qsort(rlist, (size_t) num_relations2, sizeof(relation_t), qsieve_compare_relation);
           qsieve_insert_relation(qs_inf, rlist, num_relations2);
        }
    }

    /* the full relations were moved to rlist */
    for (i = num_fulls; i < num_relations; i++)
        qsieve_relation_clear(rel_list + i);
    flint_free(rel_list);

    for (i = 0; i < rlist_length; i++)
        qsieve_relation_clear(rlist + i);
    flint_free(rlist);

    flint_free(verts);
    flint_free(eu);
    flint_free(ev);
    flint_free(alive);
    flint_free(tree);
    flint_free(deg);
    flint_free(off);
    flint_free(adj);
    flint_free(queue);
    flint_free(par_v);
    flint_free(par_e);
    flint_free(depth);

    return done;
}
