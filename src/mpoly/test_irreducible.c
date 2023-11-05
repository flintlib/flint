/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "ulong_extras.h"
#include "long_extras.h"
#include "mpoly.h"

typedef struct {
    slong x;
    slong y;
} point2d;

static int _point2d_cmp(slong x1, slong y1, slong x2, slong y2)
{
    if (x1 < x2)
        return -1;
    if (x1 > x2)
        return 1;
    if (y1 < y2)
        return -1;
    if (y1 > y2)
        return 1;
    return 0;
}

static int point2d_cmp(const point2d * p, const point2d * q)
{
    return _point2d_cmp(p->x, p->y, q->x, q->y);
}

/*
    Standing on O and looking at A, is B strictly to the left?
    i.e. (A.y - O.y) * (B.x - O.x) - (A.x - O.x) * (B.y - O.y) < 0
*/
static int _is_ccw(point2d O, point2d A, point2d B)
{
    return z_mat22_det_is_negative(A.y - O.y, A.x - O.x,
                                   B.y - O.y, B.x - O.x);
}


typedef struct {
    point2d * points;
    slong length;
    slong alloc;
} point2d_set_struct;

typedef point2d_set_struct point2d_set_t[1];

#ifdef FLINT_WANT_ASSERT
static int point2d_set_is_canonical(const point2d_set_t A)
{
    slong i;
    for (i = 1; i < A->length; i++)
    {
        if (point2d_cmp(A->points + i, A->points + i - 1) <= 0)
            return 0;
    }

    return 1;
}
#endif

static void point2d_set_init(point2d_set_t A)
{
    A->points = NULL;
    A->length = 0;
    A->alloc = 0;
}

static void point2d_set_clear(point2d_set_t A)
{
    if (A->alloc > 0)
        flint_free(A->points);
}

static void point2d_set_swap(point2d_set_t A, point2d_set_t B)
{
    point2d_set_struct t = *B;
    *B = *A;
    *A = t;
}

static void point2d_set_fit_length(point2d_set_t A, slong len)
{
    if (len > A->alloc)
    {
        len = FLINT_MAX(len, 2*A->alloc);
        A->points = (point2d *) flint_realloc(A->points, len*sizeof(point2d));
        A->alloc = len;
    }
}

static void point2d_set_sort(point2d_set_t A)
{
    qsort(A->points, A->length, sizeof(point2d),
                                 (int(*)(const void*,const void*))point2d_cmp);
}

/*
    P is a sorted array of nP distinct points
    compute the points on a ccw traversal of the convex hull via
        Andrew's monotone chain convex hull algorithm
*/
static slong convex_hull_ccw(slong * idxs, const point2d * P, slong nP)
{
    slong i, j, k = 0;

    if (nP < 3)
    {
        for (i = 0; i < nP; i++)
            idxs[i] = i;

        return nP;
    }

    for (i = 0; i < nP; i++)
    {
        while (k >= 2 && !_is_ccw(P[idxs[k - 2]], P[idxs[k - 1]], P[i]))
            k--;
        idxs[k++] = i;
    }

    for (i = nP - 1, j = k + 1; i > 0; i--)
    {
        while (k >= j && !_is_ccw(P[idxs[k - 2]], P[idxs[k - 1]], P[i - 1]))
            k--;
        idxs[k++] = i - 1;
    }

    return k - 1;
}

static int _is_in_polygon(
    point2d * V,    /* ccw polygon is V[0] -> ... -> V[nV-1] -> V[0] */
    slong nV,
    point2d p)
{
    slong i, a, b, c;
#ifdef FLINT_WANT_ASSERT
    int check;

    i = nV - 1;
    check = !_is_ccw(V[0], V[i], p);
    for (i = nV - 2; i >= 0; i--)
        check &= !_is_ccw(V[i + 1], V[i], p);
#endif

again:

    FLINT_ASSERT(nV >= 3);

    if (nV < 8)
    {
        i = nV - 1;
        if (_is_ccw(V[0], V[i], p))
        {
            FLINT_ASSERT(check == 0);
            return 0;
        }
        for (i = nV - 2; i >= 0; i--)
        {
            if (_is_ccw(V[i + 1], V[i], p))
            {
                FLINT_ASSERT(check == 0);
                return 0;
            }
        }

        FLINT_ASSERT(check == 1);
        return 1;
    }

    a = nV/4;
    b = nV/2;
    c = nV - nV/4;

    if (_is_ccw(V[a], V[0], p))
    {
        V += 0;
        nV = 1 + a;
        goto again;
    }

    if (_is_ccw(V[b], V[a], p))
    {
        V += a;
        nV = 1 + b - a;
        goto again;
    }

    if (_is_ccw(V[c], V[b], p))
    {
        V += b;
        nV = 1 + c - b;
        goto again;
    }

    if (!_is_ccw(V[0], V[c], p))
    {
        FLINT_ASSERT(check == 1);
        return 1;
    }

    if (!_is_ccw(V[nV - 1], V[c], p))
    {
        FLINT_ASSERT(check == !_is_ccw(V[0], V[nV - 1], p));
        return !_is_ccw(V[0], V[nV - 1], p);
    }

    V += c;
    nV = nV - c;

    if (nV >= 3)
        goto again;

    FLINT_ASSERT(nV == 2);
    FLINT_ASSERT(check == 0);
    return 0;
}


/* T = A union ((B + s) intersect V) */
static void point2d_set_merge_shift(
    point2d_set_t T,
    const point2d_set_t A,
    const point2d_set_t B,
    slong sx, slong sy,
    point2d * V, slong nV)
{
    point2d * Apoints, * Bpoints, * Tpoints;
    slong Alen = A->length;
    slong Blen = B->length;
    slong i, j, k;

    FLINT_ASSERT(T != A);
    FLINT_ASSERT(T != B);

    point2d_set_fit_length(T, Alen + Blen);

    Apoints = A->points;
    Bpoints = B->points;
    Tpoints = T->points;

    i = j = k = 0;

    while (i < Alen && j < Blen)
    {
        slong Bsx = Bpoints[j].x + sx;
        slong Bsy = Bpoints[j].y + sy;
        int cmp = _point2d_cmp(Apoints[i].x, Apoints[i].y, Bsx, Bsy);

        if (cmp < 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            j++;
            k++;
        }
        else
        {
            Tpoints[k].x = Bsx;
            Tpoints[k].y = Bsy;
            j++;
            k += _is_in_polygon(V, nV, Tpoints[k]);
        }

    }

    while (i < Alen)
    {
        Tpoints[k] = Apoints[i];
        i++;
        k++;
    }

    while (j < Blen)
    {
        Tpoints[k].x = Bpoints[j].x + sx;
        Tpoints[k].y = Bpoints[j].y + sy;
        j++;
        k += _is_in_polygon(V, nV, Tpoints[k]);
    }

    T->length = k;

    FLINT_ASSERT(point2d_set_is_canonical(T));
}

/* T = A union B */
static void point2d_set_merge(
    point2d_set_t T,
    const point2d_set_t A,
    const point2d_set_t B)
{
    point2d * Apoints, * Bpoints, * Tpoints;
    slong Alen = A->length;
    slong Blen = B->length;
    slong i, j, k;

    FLINT_ASSERT(T != A);
    FLINT_ASSERT(T != B);

    point2d_set_fit_length(T, Alen + Blen);

    Apoints = A->points;
    Bpoints = B->points;
    Tpoints = T->points;

    i = j = k = 0;

    while (i < Alen && j < Blen)
    {
        int cmp = _point2d_cmp(Apoints[i].x, Apoints[i].y,
                               Bpoints[j].x, Bpoints[j].y);
        if (cmp < 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            j++;
            k++;
        }
        else
        {
            Tpoints[k] = Bpoints[j];
            j++;
            k += 1;
        }
    }

    while (i < Alen)
    {
        Tpoints[k] = Apoints[i];
        i++;
        k++;
    }

    while (j < Blen)
    {
        Tpoints[k] = Bpoints[j];
        j++;
        k += 1;
    }

    T->length = k;

    FLINT_ASSERT(point2d_set_is_canonical(T));
}



#ifdef FLINT_WANT_ASSERT
static int point2d_set_contains(const point2d_set_t A, slong x, slong y)
{
    slong lo = 0;
    slong mid;
    slong hi = A->length;
    const point2d * Apoints = A->points;
    int cmp;

again:

    if (hi - lo < 8)
    {
        for ( ; lo < hi; lo++)
        {
            if (Apoints[lo].x == x && Apoints[lo].y == y)
                return 1;
        }
        return 0;
    }

    mid = lo + (hi - lo)/2;

    cmp = _point2d_cmp(Apoints[mid].x, Apoints[mid].y, x, y);

    if (cmp == 0)
        return 1;

    if (cmp < 0)
        lo = mid;
    else
        hi = mid;

    goto again;
}
#endif

/* is A intersect B empty? */
int point2d_set_disjoint(
    const point2d_set_t A,
    const point2d_set_t B)
{
    point2d * Apoints = A->points;
    point2d * Bpoints = B->points;
    slong Alen = A->length;
    slong Blen = B->length;
    slong lo, mid, hi;
    int cmp;
#ifdef FLINT_WANT_ASSERT
    int check = 1;

    for (lo = 0; lo < Blen; lo++)
        check &= !point2d_set_contains(A, Bpoints[lo].x, Bpoints[lo].y);

#endif

again:

    if (Alen < 1 || Blen < 1)
    {
        FLINT_ASSERT(check == 1);
        return 1;
    }

    if (Alen < Blen)
    {
        FLINT_SWAP(slong, Alen, Blen);
        {
            point2d * T = Apoints;
            Apoints = Bpoints;
            Bpoints = T;
        }
    }

    cmp = point2d_cmp(&Bpoints[0], &Apoints[0]);

    if (cmp == 0)
    {
        FLINT_ASSERT(check == 0);
        return 0;
    }

    if (cmp < 0)
    {
        Bpoints += 1;
        Blen -= 1;
        goto again;
    }

    /*
        throw out everything from A that is < B[0]
        if A contains B[0], return 0
    */

    lo = 0;
    hi = Alen - 1;

    cmp = point2d_cmp(&Bpoints[0], &Apoints[hi]);
    if (cmp >= 0)
    {
        FLINT_ASSERT(cmp == check);
        return cmp;
    }

search:

    FLINT_ASSERT(point2d_cmp(&Apoints[lo], &Bpoints[0]) < 0);
    FLINT_ASSERT(point2d_cmp(&Bpoints[0], &Apoints[hi]) < 0);

    if (hi - lo < 8)
    {
        for (lo++ ; lo < hi; lo++)
        {
            cmp = point2d_cmp(&Bpoints[0], &Apoints[lo]);

            if (cmp == 0)
            {
                FLINT_ASSERT(check == 0);
                return 0;
            }

            if (cmp < 0)
                break;
        }

        Apoints += lo;
        Alen -= lo;
        Bpoints += 1;
        Blen -= 1;
        goto again;
    }

    mid = lo + (hi - lo)/2;

    cmp = point2d_cmp(&Apoints[mid], &Bpoints[0]);

    if (cmp == 0)
    {
        FLINT_ASSERT(check == 0);
        return 0;
    }

    if (cmp < 0)
        lo = mid;
    else
        hi = mid;

    goto search;
}


/*
    ccw polygon is V[0] -> ... -> V[nV-1] -> V[0]
    |verts coordinates| < 2^(FLINT_BITS - 3)
*/
static int convex_hull_is_indecomposable(
    point2d * V,
    slong nV,
    ulong bound,
    point2d_set_t Ai,   /* tmp storage */
    point2d_set_t Aim1,
    point2d_set_t T,
    point2d_set_t S,
    point2d * E,        /* tmp of length nV */
    slong * Egcd)       /* tmp of length nV */
{
    slong i, j, k, g, prevx, prevy;
    ulong prod_hi, prod_lo;

    FLINT_ASSERT(nV >= 3);

    if (nV == 3)
    {
        ulong g = FLINT_ABS(V[2].x - V[0].x);
        g = n_gcd(g, FLINT_ABS(V[2].y - V[0].y));
        g = n_gcd(g, FLINT_ABS(V[1].x - V[0].x));
        g = n_gcd(g, FLINT_ABS(V[1].y - V[0].y));
        return g == 1;
    }

    /*
        since |V[i]| < 2^(FLINT_BITS - 3), all |E[i]| < 2^(SMALL_FMPZ_BITCOUNT_MAX) and
        no addition V[i] + E[j] will overflow
    */
    prevx = V[0].x;
    prevy = V[0].y;
    g = 0;
    prod_hi = 0; prod_lo = 1;
    for (i = nV - 1; i >= 0; i--)
    {
        E[i].x = prevx - V[i].x;
        E[i].y = prevy - V[i].y;
        prevx = V[i].x;
        prevy = V[i].y;
        Egcd[i] = n_gcd(FLINT_ABS(E[i].x), FLINT_ABS(E[i].y));
        E[i].x /= Egcd[i];
        E[i].y /= Egcd[i];
        g = n_gcd(g, Egcd[i]);
        umul_ppmm(prod_hi, prod_lo, prod_lo, Egcd[i]);
        if (prod_hi != 0)
            return -1;
    }

    if (g > 1)
        return 0;

    if (prod_lo > bound)
        return -1;

    /* S = {V[0] + j*E[nV-1]}_j */
    point2d_set_fit_length(S, Egcd[nV - 1]);
    for (j = 0; j < Egcd[nV - 1]; j++)
    {
        S->points[j].x = V[0].x - j*E[nV - 1].x;
        S->points[j].y = V[0].y - j*E[nV - 1].y;
    }
    S->length = Egcd[nV - 1];
    point2d_set_sort(S);

    /* A_{i-1} is empty */
    Aim1->length = 0;

    for (i = 0; i < nV - 1; i++)
    {
        point2d_set_fit_length(Ai, Egcd[i]);
        k = 0;
        for (j = 1; j <= Egcd[i]; j++)
        {
            Ai->points[k].x = V[0].x + j*E[i].x;
            Ai->points[k].y = V[0].y + j*E[i].y;
            if (!_is_in_polygon(V, nV, Ai->points[k]))
                break;
            k++;
        }
        Ai->length = k;
        point2d_set_sort(Ai);

        if (Aim1->length > 0)
        {
            point2d_set_merge(T, Ai, Aim1);
            point2d_set_swap(Ai, T);
            for (j = 1; j <= Egcd[i]; j++)
            {
                point2d_set_merge_shift(T, Ai, Aim1, j*E[i].x, j*E[i].y, V, nV);
                point2d_set_swap(Ai, T);
                if (!point2d_set_disjoint(Ai, S))
                    return 0;
            }
        }
        else
        {
            if (!point2d_set_disjoint(Ai, S))
                return 0;
        }

        point2d_set_swap(Aim1, Ai);
    }

    return 1;
}

static void z_rand_vec_primitive(
    slong * v, slong len,
    flint_rand_t state,
    mp_limb_t bound)
{
    slong i, g;

again:

    g = 0;
    for (i = 0; i < len; i++)
    {
        v[i] = z_randint(state, bound);
        g = n_gcd(g, FLINT_ABS(v[i]));
    }

    if (g == 0)
        goto again;

    if (g == 1)
        return;

    for (i = 0; i < len; i++)
        v[i] /= g;
}


static int _test_indecomposable2(
    slong * a,
    slong * b,
    slong n)
{
    slong i;
    ulong g = 0;
    for (i = 0; i < n; i++)
        g = n_gcd(g, FLINT_ABS(a[i] - b[i]));
    return g == 1;
}


static int _test_colinear(
    slong * a,
    slong * b,
    slong * c,
    slong n)
{
    slong i;
    int result;
    fmpz_t tn, td, sn, sd, g;

    fmpz_init(tn);
    fmpz_init(td);
    fmpz_init(sn);
    fmpz_init(sd);
    fmpz_init(g);

    for (i = 0; i < n; i++)
    {
        fmpz_set_si(sn, a[i]);
        fmpz_sub_si(sn, sn, c[i]);
        fmpz_set_si(sd, a[i]);
        fmpz_sub_si(sd, sd, b[i]);
        fmpz_gcd(g, sn, sd);
        if (fmpz_is_zero(g))
            continue;

        if (fmpz_sgn(sd) < 0)
            fmpz_neg(g, g);
        fmpz_divexact(sn, sn, g);
        fmpz_divexact(sd, sd, g);

        if (fmpz_is_zero(td))
        {
            fmpz_swap(td, sd);
            fmpz_swap(tn, sn);
        }
        else if (!fmpz_equal(sd, td) || !fmpz_equal(sn, tn))
        {
            result = 4;
            goto cleanup;
        }
    }

    if (fmpz_is_zero(td))
        result = 0;
    else if (fmpz_sgn(tn) < 0)
        result = 1;
    else if (fmpz_cmp(tn, td) > 0)
        result = 2;
    else
        result = 3;

cleanup:

    fmpz_clear(tn);
    fmpz_clear(td);
    fmpz_clear(sn);
    fmpz_clear(sd);
    fmpz_clear(g);

    return result;
}

static int _test_indecomposable3(
    slong * a,
    slong * b,
    slong * c,
    slong n)
{
    slong i;
    ulong g;

    switch (_test_colinear(a, b, c, n))
    {
        case 0:
            return 0;
        case 1:
            return _test_indecomposable2(c, b, n);
        case 2:
            return _test_indecomposable2(a, c, n);
        case 3:
            return _test_indecomposable2(a, b, n);
        default:
            break;
    }

    g = 0;
    for (i = 0; i < n; i++)
    {
        g = n_gcd(g, FLINT_ABS(a[i] - b[i]));
        g = n_gcd(g, FLINT_ABS(a[i] - c[i]));
    }
    return g == 1;
}

/*
    Fast Absolute Irreducibility Testing via Newton Polytopes
    Shuhong Gao and Alan G.B. Lauderz

    Aexps is the Alen x nvars exponent matrix with the entry A[i,j]
    at Aexps[i*stride + j]
*/
int _mpoly_test_irreducible(
    slong * Aexps, slong stride, slong Alen,
    slong nvars,
    flint_rand_t state,
    slong tries_left)   /* what they call the "projection bound" */
{
    int success;
    slong i, j, newlen, hull_len;
    ulong matrix_bound = 2;
    ulong memory_bound = 1000;
    ulong max_memory_bound = UWORD(1) << (20 + FLINT_BITS/8);
    slong * rowx;
    slong * rowy;
    point2d * points;
    slong * dups;
    slong * hull_idxs;
    slong hull_idxs_alloc;
    point2d * hull_points;
    slong hull_points_alloc;
    point2d_set_t T1, T2, T3, T4;

    if (Alen < 2 || nvars < 2)
        return 0;

    if (Alen == 2)
        return _test_indecomposable2(Aexps + 0*stride, Aexps + 1*stride, nvars);

    if (Alen == 3)
        return _test_indecomposable3(Aexps + 0*stride, Aexps + 1*stride,
                                                      Aexps + 2*stride, nvars);
    if (tries_left <= 0)
        return 0;

    rowx = FLINT_ARRAY_ALLOC(2*nvars, slong);
    rowy = rowx + nvars;
    points = FLINT_ARRAY_ALLOC(Alen, point2d);
    dups = FLINT_ARRAY_ALLOC(Alen, slong);
    hull_idxs_alloc = 10;
    hull_idxs = FLINT_ARRAY_ALLOC(hull_idxs_alloc, slong);
    hull_points_alloc = 4;
    hull_points = FLINT_ARRAY_ALLOC(hull_points_alloc, point2d);

    point2d_set_init(T1);
    point2d_set_init(T2);
    point2d_set_init(T3);
    point2d_set_init(T4);

again:

    if (--tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    memory_bound = FLINT_MIN(max_memory_bound, memory_bound/8*9);
    matrix_bound += 1;

    if (nvars == 2)
    {
        tries_left = 0;
        memory_bound = max_memory_bound;

        for (i = 0; i < Alen; i++)
        {
            slong lox = Aexps[i*stride + 0];
            slong loy = Aexps[i*stride + 1];

            if (FLINT_MIN(lox, loy) <= -(WORD(1) << (FLINT_BITS - 3)) ||
                FLINT_MAX(lox, loy) >= WORD(1) << (FLINT_BITS - 3))
            {
                success = 0;
                goto cleanup;
            }

            points[i].x = lox;
            points[i].y = loy;
        }
    }
    else
    {
        z_rand_vec_primitive(rowx, nvars, state, matrix_bound);
        z_rand_vec_primitive(rowy, nvars, state, matrix_bound);

        for (i = 0; i < Alen; i++)
        {
            ulong x2, x1, x0, y2, y1, y0, p2, p1, p0;
            x2 = x1 = x0 = y2 = y1 = y0 = 0;
            for (j = 0; j < nvars; j++)
            {
                smul_ppmm(p1, p0, Aexps[i*stride + j], rowx[j]);
                p2 = FLINT_SIGN_EXT(p1);
                add_sssaaaaaa(x2, x1, x0, x2, x1, x0, p2, p1, p0);

                smul_ppmm(p1, p0, Aexps[i*stride + j], rowy[j]);
                p2 = FLINT_SIGN_EXT(p1);
                add_sssaaaaaa(y2, y1, y0, y2, y1, y0, p2, p1, p0);
            }

            if (x2 != FLINT_SIGN_EXT(x0) || x1 != FLINT_SIGN_EXT(x0) ||
                y2 != FLINT_SIGN_EXT(y0) || y1 != FLINT_SIGN_EXT(y0))
            {
                goto again;
            }

            points[i].x = x0;
            points[i].y = y0;

            if (points[i].x <= -(WORD(1) << (FLINT_BITS - 3)) ||
                points[i].y <= -(WORD(1) << (FLINT_BITS - 3)) ||
                points[i].x >= WORD(1) << (FLINT_BITS - 3) ||
                points[i].y >= WORD(1) << (FLINT_BITS - 3))
            {
                goto again;
            }
        }
    }

    qsort(points, Alen, sizeof(point2d),
                                 (int(*)(const void*,const void*))point2d_cmp);

    /* delete duplicates and track which are duplicated */
    dups[0] = 0;
    newlen = 1;
    for (i = 1; i < Alen; i++)
    {
        if (point2d_cmp(&points[newlen - 1], &points[i]) == 0)
        {
            dups[newlen - 1] = 1;
        }
        else
        {
            dups[newlen] = 0;
            points[newlen] = points[i];
            newlen++;
        }
    }

    /* find indices of convex hull */
    if (newlen + 1 > hull_idxs_alloc)
    {
        hull_idxs_alloc = newlen + 1;
        hull_idxs = (slong *) flint_realloc(hull_idxs,
                                                hull_idxs_alloc*sizeof(slong));
    }
    hull_len = convex_hull_ccw(hull_idxs, points, newlen);
    if (hull_len < 3)
        goto again;

    /* ensure no duplicates on hull */
    if (hull_len + 1 > hull_points_alloc)
    {
        hull_points_alloc = hull_len + 1;
        hull_points = (point2d *) flint_realloc(hull_points,
                                            hull_points_alloc*sizeof(point2d));
    }
    for (i = 0; i < hull_len; i++)
    {
        hull_points[i] = points[hull_idxs[i]];
        if (dups[hull_idxs[i]] != 0)
            goto again;
    }

    /* check indecomposability with a bound on the memory usage */
    success = convex_hull_is_indecomposable(hull_points, hull_len,
                                   memory_bound, T1, T2, T3, T4, points, dups);
    if (success < 1)
    {
        if (success < 0)
            memory_bound = FLINT_MIN(max_memory_bound, memory_bound/8*9);
        goto again;
    }

    success = 1;

cleanup:

    point2d_set_clear(T1);
    point2d_set_clear(T2);
    point2d_set_clear(T3);
    point2d_set_clear(T4);

    flint_free(rowx);
    flint_free(points);
    flint_free(dups);
    flint_free(hull_idxs);
    flint_free(hull_points);

    return success;
}


int mpoly_test_irreducible(
    ulong * Aexps, flint_bitcnt_t Abits, slong Alen,
    const mpoly_ctx_t ctx)
{
    int result, overflowed;
    slong n = ctx->nvars;
    slong i, j, N, sum_deg, tries;
    slong * uexps, * max_exps;
    flint_rand_t state;

    if (Abits > FLINT_BITS || Alen < 2)
        return 0;

    flint_randinit(state);

    uexps = FLINT_ARRAY_ALLOC(n*Alen, slong);
    max_exps = FLINT_ARRAY_ALLOC(n, slong);

    for (j = 0; j < n; j++)
        max_exps[j] = 0;

    N = mpoly_words_per_exp(Abits, ctx);

    for (i = 0; i < Alen; i++)
    {
        mpoly_get_monomial_ui_sp((ulong *)uexps + n*i, Aexps + N*i, Abits, ctx);
        for (j = 0; j < n; j++)
            max_exps[j] = FLINT_MAX(max_exps[j], uexps[n*i + j]);
    }

    sum_deg = 1;
    overflowed = 0;
    for (j = 0; j < n; j++)
    {
        if (z_add_checked(&sum_deg, sum_deg, max_exps[j]))
        {
            overflowed = 1;
            break;
        }
    }

    tries = 12;
    if (!overflowed)
        tries -= Alen/sum_deg/2;

    result = _mpoly_test_irreducible(uexps, n, Alen, n, state, tries);

    flint_randclear(state);
    flint_free(max_exps);
    flint_free(uexps);

    return result;
}

