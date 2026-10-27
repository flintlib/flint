/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2021, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "mpn_extras.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"

/*
    The content of a vector (see doc/src/mpn_extras/gcd_vec.md).

    After the two end entries, the entries are taken in blocks alternately
    from the front and the back of the vector (the leading and trailing
    coefficients of polynomials being the most likely to differ from the
    rest in their factors), with sizes doubling from 4 so that little is
    wasted when the content is found to be 1 early. The small entries of a
    block are passed in a buffer to flint_mpn_gcd_vec_1, and its large ones
    by pointer to flint_mpn_gcd_vec, which is cheap once the running gcd
    has a single limb.

    A gcd of large entries is expensive, and can be avoided altogether if
    small entries (or shorter large entries) come later. So when the first
    large entry is long compared with the length of the vector, the large
    entries are deferred to a second pass: the first pass does the small
    entries and finds the two shortest large entries, from which the second
    pass starts if there are no small entries. Otherwise (entries short
    compared with the length) each large entry is done in its block.
*/

#define CONTENT_BLOCK_SMALL 256
#define CONTENT_BLOCK_BIG 64

/* defer the large entries if their length n satisfies
   n^2 >= CONTENT_DEFER_RATIO * (number of entries left) */
#define CONTENT_DEFER_RATIO 1

#define BLOCK_NEXT(s, e, lo, hi, b, side, bmax)         \
    do {                                                \
        if ((side) == 0)                                \
        {                                               \
            (s) = (lo);                                 \
            (e) = (lo) = FLINT_MIN((lo) + (b), (hi));   \
        }                                               \
        else                                            \
        {                                               \
            (e) = (hi);                                 \
            (s) = (hi) = FLINT_MAX((hi) - (b), (lo));   \
            (b) = FLINT_MIN(2 * (b), (bmax));           \
        }                                               \
        (side) ^= 1;                                    \
    } while (0)

#define CONTENT_STACK_LIMBS 16

/* The running gcd: g if gp == NULL (g = 0 if no entries yet), else
   {gp, gn} with gn >= 2, where gp is stack (if short enough) or heap. */
typedef struct
{
    ulong g;
    mp_ptr gp;
    mp_size_t gn;
    mp_limb_t stack[CONTENT_STACK_LIMBS];
}
content_state_t;

static void
_content_free(content_state_t * S)
{
    if (S->gp != NULL && S->gp != S->stack)
        flint_free(S->gp);
    S->gp = NULL;
}

/* replace the running gcd by its gcd with the given large entries;
   returns 1 if the result is 1 */
static int
_content_big(content_state_t * S, mp_srcptr * xp, mp_size_t * xn, slong n)
{
    mp_limb_t g1;
    slong i;

    if (n == 0)
        return 0;

    if (S->gp == NULL && S->g != 0)
    {
        g1 = S->g;
        flint_mpn_gcd_vec(&g1, 1, xp, xn, n);
        S->g = g1;
        return g1 == 1;
    }

    if (S->gp == NULL)
    {
        /* the gcd fits in the shortest of these entries */
        mp_size_t nmin = xn[0];
        for (i = 1; i < n; i++)
            nmin = FLINT_MIN(nmin, xn[i]);
        S->gp = (nmin <= CONTENT_STACK_LIMBS) ? S->stack : flint_malloc(nmin * sizeof(mp_limb_t));
        S->gn = 0;
    }

    S->gn = flint_mpn_gcd_vec(S->gp, S->gn, xp, xn, n);

    if (S->gn == 1)
    {
        S->g = S->gp[0];
        _content_free(S);
        return S->g == 1;
    }

    return 0;
}

/* replace the running gcd by its gcd with the given single limbs */
static int
_content_small(content_state_t * S, mp_srcptr x, slong n)
{
    if (S->gp == NULL)
    {
        S->g = flint_mpn_gcd_vec_1(x, n, S->g);
    }
    else
    {
        mp_limb_t g1 = flint_mpn_gcd_vec_1(x, n, 0);
        mp_srcptr gp = S->gp;

        if (g1 == 0)
            return 0;

        flint_mpn_gcd_vec(&g1, 1, &gp, &S->gn, 1);
        _content_free(S);
        S->g = g1;
    }

    return S->g == 1;
}

#define BIG_PTR(xp, xn, c)                      \
    do {                                        \
        mpz_ptr __z = COEFF_TO_PTR(c);          \
        (xp) = __z->_mp_d;                      \
        (xn) = FLINT_ABS(__z->_mp_size);        \
    } while (0)

/* Short vectors: the block structure costs more than it saves, so all
   small entries go to a single flint_mpn_gcd_vec_1 call (which takes
   them from both ends) and all large entries to a single
   flint_mpn_gcd_vec call (which starts from the shortest one). */
#define CONTENT_SHORT 16

static void
_fmpz_vec_content_short(fmpz_t res, const fmpz * vec, slong len, ulong g0, const fmpz * inp)
{
    mp_limb_t buf[CONTENT_SHORT];
    mp_srcptr xp[CONTENT_SHORT + 1];
    mp_size_t xn[CONTENT_SHORT + 1];
    mp_limb_t g1;
    slong i, j, nb;
    mp_size_t gn, nmin;

    g1 = g0;
    nb = 0;
    if (inp != NULL)
    {
        BIG_PTR(xp[0], xn[0], *inp);
        nb = 1;
    }

    /* the entries alternately from both ends */
    for (i = j = 0; i < len; i++)
    {
        fmpz c = vec[(i & 1) ? len - 1 - i / 2 : i / 2];
        if (!COEFF_IS_MPZ(c))
            buf[j++] = FLINT_ABS(c);
        else
        {
            BIG_PTR(xp[nb], xn[nb], c);
            nb++;
        }
    }

    if (j != 0)
    {
        g1 = flint_mpn_gcd_vec_1(buf, j, g1);
        if (g1 == 1 || nb == 0)
        {
            fmpz_set_ui(res, g1);
            return;
        }
    }

    if (g1 != 0)
    {
        flint_mpn_gcd_vec(&g1, 1, xp, xn, nb);
        fmpz_set_ui(res, g1);
        return;
    }

    nmin = xn[0];
    for (i = 1; i < nb; i++)
        nmin = FLINT_MIN(nmin, xn[i]);

    {
        mp_ptr gp;
        mpz_ptr z;
        TMP_INIT;
        TMP_START;
        gp = TMP_ALLOC(nmin * sizeof(mp_limb_t));
        gn = flint_mpn_gcd_vec(gp, 0, xp, xn, nb);
        if (gn == 1)
        {
            fmpz_set_ui(res, gp[0]);
        }
        else
        {
            z = _fmpz_promote(res);
            flint_mpn_copyi(FLINT_MPZ_REALLOC(z, gn), gp, gn);
            z->_mp_size = gn;
        }
        TMP_END;
    }
}

static void
_fmpz_vec_content_long(fmpz_t res, const fmpz * vec, slong len, ulong g0, const fmpz * inp)
{
    mp_limb_t buf[CONTENT_BLOCK_SMALL];
    mp_srcptr xp[CONTENT_BLOCK_BIG + 1];
    mp_size_t xn[CONTENT_BLOCK_BIG + 1];
    content_state_t S[1];
    slong i, j, nb, s, e, lo, hi, b, fl, ndefer = 0, imin = -1, imin2 = -1;
    slong nmin = WORD_MAX, nmin2 = WORD_MAX;   /* mp_size_t may be 32 bits */
    mp_size_t n;
    int side, defer = -1;
    /* set on first use (tg = 0 differs from any g > 1) */
    ulong tg, th = 1, thinv = 1, tlim = UWORD_MAX, tmask = 0;
    unsigned int tk;
    int prev_big = 0;

    S->g = g0;
    S->gp = NULL;
    S->gn = 0;
    nb = 0;

    if (inp != NULL)
    {
        BIG_PTR(xp[0], xn[0], *inp);
        nb = 1;
    }

    lo = 0;
    hi = len;

    /* first pass; the large entries done in it are passed in groups of
       at least fl (doubling up to CONTENT_BLOCK_BIG) */
    b = 4;
    side = 0;
    fl = 2;
    tg = 0;
    while (lo < hi)
    {
        BLOCK_NEXT(s, e, lo, hi, b, side, CONTENT_BLOCK_SMALL);

        /* with a single-limb running gcd 2^k h, a block of small entries
           all divisible by it is common and is detected directly on the
           fmpz words: |c| has its low k bits zero and |c| h^{-1} mod B
           <= (B - 1) / h, and c is not a pointer; it is only tried after a
           block without large entries */
        if (S->gp == NULL && S->g > 1 && !prev_big)
        {
            ulong bad = 0, a, sg;
            fmpz c;

            if (S->g != tg)
            {
                tg = S->g;
                tk = flint_ctz(tg);
                tmask = (UWORD(1) << tk) - 1;
                th = tg >> tk;
                thinv = n_binvert(th);
                tlim = UWORD_MAX / th;
            }

            for (i = s; i < e; i++)
            {
                c = vec[i];
                sg = FLINT_SIGN_EXT(c);
                a = (c ^ sg) - sg;
                bad |= (a & tmask) | (ulong) (a * thinv > tlim) | (ulong) COEFF_IS_MPZ(c);
            }

            if (bad == 0)
                continue;
        }

        for (i = s, j = 0; i < e; i++)
        {
            fmpz c = vec[i];

            if (!COEFF_IS_MPZ(c))
            {
                buf[j++] = FLINT_ABS(c);
            }
            else
            {
                BIG_PTR(xp[nb], xn[nb], c);

                if (FLINT_UNLIKELY(defer == -1))
                {
                    /* the cost of a gcd against that of the rest of
                       the pass, if it turns out to be unnecessary */
                    double r = (double) xn[nb] * (double) xn[nb];
                    defer = (r >= CONTENT_DEFER_RATIO * (double) (hi - lo + e - i));
                }

                if (defer)
                {
                    ndefer++;
                    if (xn[nb] < nmin2)
                    {
                        if (xn[nb] < nmin)
                        {
                            nmin2 = nmin;
                            imin2 = imin;
                            nmin = xn[nb];
                            imin = i;
                        }
                        else
                        {
                            nmin2 = xn[nb];
                            imin2 = i;
                        }
                    }
                }
                else
                {
                    nb++;
                    if (nb == CONTENT_BLOCK_BIG)
                    {
                        if (_content_big(S, xp, xn, nb))
                            goto one;
                        nb = 0;
                    }
                }
            }
        }

        prev_big = (j != e - s);

        if (_content_small(S, buf, j))
            goto one;

        if (nb >= fl && defer == 0)
        {
            if (_content_big(S, xp, xn, nb))
                goto one;
            nb = 0;
            fl = FLINT_MIN(2 * fl, CONTENT_BLOCK_BIG);
        }
    }

    if (nb != 0)
    {
        if (_content_big(S, xp, xn, nb))
            goto one;
        nb = 0;
    }

    /* second pass over the deferred large entries, starting from the two
       shortest ones if the running gcd is not a single limb */
    if (ndefer != 0)
    {
        if (S->gp == NULL && S->g == 0)
        {
            BIG_PTR(xp[0], xn[0], vec[imin]);
            nb = 1;
            if (imin2 != -1)
            {
                BIG_PTR(xp[1], xn[1], vec[imin2]);
                nb = 2;
            }
            if (_content_big(S, xp, xn, nb))
                goto one;
            nb = 0;
        }
        else
        {
            imin = imin2 = -1;
        }

        lo = 0;
        hi = len;
        b = 1;
        side = 0;
        fl = 4;
        while (lo < hi)
        {
            BLOCK_NEXT(s, e, lo, hi, b, side, CONTENT_BLOCK_BIG);

            for (i = s; i < e; i++)
            {
                if (COEFF_IS_MPZ(vec[i]) && i != imin && i != imin2)
                {
                    BIG_PTR(xp[nb], xn[nb], vec[i]);
                    nb++;
                    if (nb == CONTENT_BLOCK_BIG)
                    {
                        if (_content_big(S, xp, xn, nb))
                            goto one;
                        nb = 0;
                    }
                }
            }

            if (nb >= fl || lo >= hi)
            {
                if (_content_big(S, xp, xn, nb))
                    goto one;
                nb = 0;
                fl = FLINT_MIN(2 * fl, CONTENT_BLOCK_BIG);
            }
        }
    }

    if (S->gp == NULL)
    {
        fmpz_set_ui(res, S->g);
    }
    else
    {
        mpz_ptr z = _fmpz_promote(res);
        n = S->gn;
        flint_mpn_copyi(FLINT_MPZ_REALLOC(z, n), S->gp, n);
        z->_mp_size = n;
        _content_free(S);
    }
    return;

one:
    _content_free(S);
    fmpz_one(res);
}

/* res = gcd(vec, init) for len >= 1, without the multithreading */
static void
_fmpz_vec_content_seq(fmpz_t res, const fmpz * vec, slong len, const fmpz_t init)
{
    ulong g0 = 0;
    const fmpz * inp = NULL;

    if (!COEFF_IS_MPZ(*init))
        g0 = FLINT_ABS(*init);
    else
        inp = init;

    if (len <= CONTENT_SHORT)
        _fmpz_vec_content_short(res, vec, len, g0, inp);
    else
        _fmpz_vec_content_long(res, vec, len, g0, inp);
}

/*
    Multithreading: once the entries at both ends (about an eighth of the
    vector, at most CONTENT_THREAD_WARMUP from each end) have been done
    without finding the content to be 1, the running gcd is probably the
    content or close to it, and what remains is mostly to check that it
    divides the other entries, which is independent work. If the other
    entries are long enough, they are split into chunks of about equal
    work whose contents (starting from the running gcd) are computed in
    parallel. The work of an entry of n limbs is estimated in units of the
    cost of reducing a limb modulo a single limb: n for a single-limb gcd,
    and max(n, (n - gn + 1) gn / 4) for a running gcd of gn limbs (the
    cost of a division). Each thread is given at least
    CONTENT_THREAD_WORK / 2 units, as starting threads costs about as much
    as reducing 10^5 limbs modulo a single limb, and the number of threads
    is at most sqrt(work / CONTENT_THREAD_START), as waking each thread
    costs about CONTENT_THREAD_START units. Before any of this, an upper
    bound for the work from the sizes of the end entries decides whether
    threads are possible at all, so that short vectors and vectors of small
    entries take exactly the single-threaded path.
*/
#define CONTENT_THREAD_MIN_LEN 64
#define CONTENT_THREAD_WARMUP 16
#define CONTENT_THREAD_WORK 600000
#define CONTENT_THREAD_START 10000

static double
_content_work(double n, double gn)
{
    if (gn <= 1)
        return n;
    return FLINT_MAX(n, (FLINT_MAX(n - gn, 0) + 1) * gn / 4);
}

/* an upper bound for the work, assuming entries of the size of the larger
   end entry and a running gcd of that size (or of a single limb if g0 is
   nonzero, as it then divides it) */
static double
_content_work_bound(const fmpz * vec, slong len, ulong g0)
{
    double s0, s1, s;

    s0 = COEFF_IS_MPZ(vec[0]) ? FLINT_ABS(COEFF_TO_PTR(vec[0])->_mp_size) : 1;
    s1 = COEFF_IS_MPZ(vec[len - 1]) ? FLINT_ABS(COEFF_TO_PTR(vec[len - 1])->_mp_size) : 1;
    s = FLINT_MAX(s0, s1);

    return (double) len * ((g0 != 0) ? s : FLINT_MAX(s, s * s / 4));
}

typedef struct
{
    const fmpz * vec;
    const slong * start;    /* chunk i is vec[start[i]:start[i + 1]] */
    const fmpz * init;
    fmpz * res;
}
content_thread_arg_t;

static void
_content_worker(slong i, void * args)
{
    content_thread_arg_t * A = args;
    slong s = A->start[i], e = A->start[i + 1];

    if (e > s)
        _fmpz_vec_content_seq(A->res + i, A->vec + s, e - s, A->init);
    else
        fmpz_abs(A->res + i, A->init);
}

/* returns 0 if threads were not used (res not set) */
static int
_fmpz_vec_content_threaded(fmpz_t res, const fmpz * vec, slong len, const fmpz_t init)
{
    slong w, i, j, k, nt, lo, hi;
    double work, acc, per, gn;
    fmpz_t g;
    fmpz * r;
    slong * start;
    content_thread_arg_t A;

    nt = flint_get_num_threads();
    if (nt <= 1 || len < CONTENT_THREAD_MIN_LEN)
        return 0;

    /* the entries at both ends, in one call (which does the small and the
       shortest large entries first), through a shallow copy of the fmpz
       words: the copy is only read and is freed without clearing */
    w = FLINT_MIN(CONTENT_THREAD_WARMUP, len / 16);
    fmpz_init(g);
    {
        fmpz sample[2 * CONTENT_THREAD_WARMUP];
        for (i = 0; i < w; i++)
        {
            sample[2 * i] = vec[i];
            sample[2 * i + 1] = vec[len - 1 - i];
        }
        _fmpz_vec_content_seq(g, sample, 2 * w, init);
    }
    if (fmpz_is_one(g))
    {
        fmpz_swap(res, g);
        fmpz_clear(g);
        return 1;
    }

    lo = w;
    hi = len - w;

#define ENTRY_LIMBS(c) (COEFF_IS_MPZ(c) ? FLINT_ABS(COEFF_TO_PTR(c)->_mp_size) : 1)

    gn = COEFF_IS_MPZ(*g) ? FLINT_ABS(COEFF_TO_PTR(*g)->_mp_size) : 1;

    /* estimate the work from the entries at both ends before scanning
       the others, as the scan would cost too much when there is not
       enough work for threads */
    work = 0;
    for (i = 0; i < w; i++)
        work += _content_work(ENTRY_LIMBS(vec[i]), gn) + _content_work(ENTRY_LIMBS(vec[len - 1 - i]), gn);
    work = work * (hi - lo) / (2 * w);

    if (work < CONTENT_THREAD_WORK / 2)
        nt = 1;

    if (nt > 1)
    {
        work = 0;
        for (i = lo; i < hi; i++)
            work += _content_work(ENTRY_LIMBS(vec[i]), gn);
        nt = FLINT_MIN(nt, (slong) (work / (CONTENT_THREAD_WORK / 2)));
        /* waking nt threads costs about nt CONTENT_THREAD_START units,
           which balances the gain at about sqrt(work / START) threads */
        nt = FLINT_MIN(nt, (slong) sqrt(work / CONTENT_THREAD_START));
    }

    nt = FLINT_MIN(nt, hi - lo);

    if (nt <= 1)
    {
        _fmpz_vec_content_seq(g, vec + lo, hi - lo, g);
        fmpz_swap(res, g);
        fmpz_clear(g);
        return 1;
    }

    /* chunks of about work / nt units */
    start = flint_malloc((nt + 1) * sizeof(slong));
    r = _fmpz_vec_init(nt);
    per = work / nt;
    start[0] = lo;
    acc = 0;
    for (i = lo, k = 1; i < hi && k < nt; i++)
    {
        acc += _content_work(ENTRY_LIMBS(vec[i]), gn);
        if (acc >= k * per)
            start[k++] = i + 1;
    }
    for ( ; k <= nt; k++)
        start[k] = hi;

    A.vec = vec;
    A.start = start;
    A.init = g;
    A.res = r;
    flint_parallel_do(_content_worker, &A, nt, nt, FLINT_PARALLEL_UNIFORM);

    /* the chunk contents divide g; take the smallest first */
    j = 0;
    for (i = 1; i < nt; i++)
        if (fmpz_cmpabs(r + i, r + j) < 0)
            j = i;
    fmpz_swap(res, r + j);
    for (i = 0; i < nt && !fmpz_is_one(res); i++)
        if (i != j)
            fmpz_gcd(res, res, r + i);

    _fmpz_vec_clear(r, nt);
    flint_free(start);
    fmpz_clear(g);
    return 1;

#undef ENTRY_LIMBS
}

/* res = gcd(vec, g0, inp) for len >= 1 nonzero end entries, g0 a single
   limb and inp NULL or a large value */
/* entries done by direct gcds after the end entries */
#define CONTENT_WARMUP 2

static void
_fmpz_vec_content_main(fmpz_t res, const fmpz * vec, slong len, ulong g0, const fmpz * inp)
{
    slong i;

    /* the small end entries first (often finding the content to be 1),
       and then a large inp against the resulting single limb, so that the
       rest is done against a single-limb gcd */
    {
        int small0 = !COEFF_IS_MPZ(vec[0]);
        int small1 = !COEFF_IS_MPZ(vec[len - 1]);

        if (small0 || small1)
        {
            mp_limb_t a, c;

            a = small0 ? FLINT_ABS(vec[0]) : FLINT_ABS(vec[len - 1]);

            /* the other end, also if it is large (one remainder pass);
               mpn_gcd_1 reduces unbalanced operands by a division first */
            if (len > 1)
            {
                if (small0 && small1)
                {
                    c = FLINT_ABS(vec[len - 1]);
                    a = mpn_gcd_1(&a, 1, c);
                }
                else
                {
                    mpz_ptr z = COEFF_TO_PTR(small0 ? vec[len - 1] : vec[0]);
                    a = mpn_gcd_1(z->_mp_d, FLINT_ABS(z->_mp_size), a);
                    small0 = small1 = 1;
                }
            }

            /* the given values last, as they often share factors with the
               entries by construction (a denominator, a previous content) */
            if (g0 != 0)
                a = mpn_gcd_1(&a, 1, g0);

            if (a != 1 && inp != NULL)
            {
                mpz_ptr z = COEFF_TO_PTR(*inp);
                a = mpn_gcd_1(z->_mp_d, FLINT_ABS(z->_mp_size), a);
                inp = NULL;
            }

            if (a == 1)
            {
                fmpz_one(res);
                return;
            }

            g0 = a;
            if (small0)
            {
                vec++;
                len--;
            }
            if (small1 && len > 0)
                len--;

            /* the next entries from both ends by direct gcds, which often
               end the computation and cost less than setting up the
               divisibility tests */
            for (i = 0; i < CONTENT_WARMUP && len > 0; i++)
            {
                fmpz e = (i & 1) ? vec[--len] : vec[0];
                if (!(i & 1))
                {
                    vec++;
                    len--;
                }
                if (!COEFF_IS_MPZ(e))
                {
                    if (e != 0)
                    {
                        c = FLINT_ABS(e);
                        a = mpn_gcd_1(&c, 1, a);
                    }
                }
                else
                {
                    mpz_ptr z = COEFF_TO_PTR(e);
                    a = mpn_gcd_1(z->_mp_d, FLINT_ABS(z->_mp_size), a);
                }
                if (a == 1)
                {
                    fmpz_one(res);
                    return;
                }
            }

            g0 = a;

            if (len <= 0)
            {
                fmpz_set_ui(res, g0);
                return;
            }
        }
    }

    if (len <= CONTENT_SHORT)
    {
        _fmpz_vec_content_short(res, vec, len, g0, inp);
    }
    else if (len >= CONTENT_THREAD_MIN_LEN && flint_get_num_threads() > 1
                && _content_work_bound(vec, len, g0) >= CONTENT_THREAD_WORK)
    {
        fmpz_t init;
        int done;

        fmpz_init(init);
        if (inp != NULL)
            fmpz_gcd_ui(init, inp, g0);   /* g0 = 0 gives |inp| */
        else
            fmpz_set_ui(init, g0);
        done = _fmpz_vec_content_threaded(res, vec, len, init);
        if (!done)
            _fmpz_vec_content_seq(res, vec, len, init);
        fmpz_clear(init);
    }
    else
    {
        _fmpz_vec_content_long(res, vec, len, g0, inp);
    }
}

void
_fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len)
{
    while (len > 0 && vec[0] == 0)
    {
        vec++;
        len--;
    }

    while (len > 0 && vec[len - 1] == 0)
        len--;

    if (len <= 2)
    {
        if (len == 0)
            fmpz_zero(res);
        else if (len == 1)
            fmpz_abs(res, vec);
        else
            fmpz_gcd(res, vec, vec + 1);
    }
    else
    {
        _fmpz_vec_content_main(res, vec, len, 0, NULL);
    }
}

void
_fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len, const fmpz_t inp)
{
    while (len > 0 && vec[0] == 0)
    {
        vec++;
        len--;
    }

    while (len > 0 && vec[len - 1] == 0)
        len--;

    if (len == 0)
        fmpz_abs(res, inp);
    else if (fmpz_is_pm1(inp))
        fmpz_one(res);
    else if (len == 1)
        fmpz_gcd(res, vec, inp);
    else if (!COEFF_IS_MPZ(*inp))
        _fmpz_vec_content_main(res, vec, len, FLINT_ABS(*inp), NULL);
    else
        _fmpz_vec_content_main(res, vec, len, 0, inp);
}
