/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "thread_support.h"
#include "mpn_extras/impl.h"

/* Default tuning parameters (see tune/profile programs). */
#ifndef FLINT_MPN_CRT_CHUNK_BITS
# define FLINT_MPN_CRT_CHUNK_BITS 2048
#endif
#ifndef FLINT_MPN_CRT_MOD_BASE_BITS
# define FLINT_MPN_CRT_MOD_BASE_BITS 2048
#endif
#ifndef FLINT_MPN_CRT_PREINV_CUTOFF
# define FLINT_MPN_CRT_PREINV_CUTOFF 16
#endif

/* parallel helper: run f over n nodes with limb size limbs, using threads
   only when the nodes are large */
static void
_crt_parallel_do(do_func_t f, void * arg, slong n, slong limbs)
{
    slong i;

    if (_crt_use_threads(n, limbs))
        flint_parallel_do(f, arg, n, -1, FLINT_PARALLEL_UNIFORM);
    else
        for (i = 0; i < n; i++)
            f(i, arg);
}


typedef struct
{
    flint_mpn_crt_struct * C;
    slong k;
    nn_ptr tmp;
}
_crt_prod_arg;

/* product of the children of node i at level k */
static void
_crt_prod_worker(slong i, void * varg)
{
    _crt_prod_arg * arg = (_crt_prod_arg *) varg;
    flint_mpn_crt_struct * C = arg->C;
    slong k = arg->k;
    slong pl = C->level_limbs[k - 1];
    slong pc = C->level_count[k - 1];
    nn_srcptr pp = C->level_prod[k - 1];
    nn_ptr v = arg->tmp + i * (2 * pl);
    slong bn = C->level_len[k - 1][2 * i];
    slong vn;

    if (2 * i + 1 < pc)
    {
        slong cn = C->level_len[k - 1][2 * i + 1];
        if (bn >= cn)
            flint_mpn_mul(v, pp + (2 * i) * pl, bn, pp + (2 * i + 1) * pl, cn);
        else
            flint_mpn_mul(v, pp + (2 * i + 1) * pl, cn, pp + (2 * i) * pl, bn);
        vn = bn + cn;
        vn -= (v[vn - 1] == 0);
    }
    else
    {
        flint_mpn_copyi(v, pp + (2 * i) * pl, bn);
        vn = bn;
    }

    C->level_len[k][i] = vn;
}

/* remainder tree step: R for the children of node i at level k */
typedef struct
{
    flint_mpn_crt_struct * C;
    slong k;
    nn_srcptr R;
    slong rl;
    nn_ptr R2;
    slong rl2;
    int bad;
}
_crt_rem_arg;

static void
_crt_rem_worker(slong i, void * varg)
{
    _crt_rem_arg * arg = (_crt_rem_arg *) varg;
    flint_mpn_crt_struct * C = arg->C;
    slong k = arg->k;
    slong pl = C->level_limbs[k - 1];
    slong pc = C->level_count[k - 1];
    nn_srcptr pp = C->level_prod[k - 1];
    nn_srcptr r = arg->R + i * arg->rl;
    slong rn = arg->rl;
    slong rl2 = arg->rl2;
    slong c;
    nn_ptr T, U, Q;

    MPN_NORM(r, rn);

    /* t: pl + 1, u: 2 pl + 1, q: rn + 1 */
    T = FLINT_ARRAY_ALLOC(4 * pl + rn + 4, ulong);
    U = T + pl + 1;
    Q = U + 2 * pl + 1;

    for (c = 2 * i; c < FLINT_MIN(2 * i + 2, pc); c++)
    {
        nn_ptr rc = arg->R2 + c * rl2;
        nn_srcptr d = pp + c * pl;
        slong dn = C->level_len[k - 1][c];
        slong o = c ^ 1;
        slong tn, un;
        nn_srcptr vo;
        slong von;

        flint_mpn_zero(rc, rl2);

        if (o >= pc)
        {
            /* single child: same product, same R */
            flint_mpn_copyi(rc, r, rn);
            continue;
        }

        /* t = r mod d */
        if (rn > dn || (rn == dn && mpn_cmp(r, d, dn) >= 0))
        {
            mpn_tdiv_qr(Q, T, 0, r, rn, d, dn);
            tn = dn;
            MPN_NORM(T, tn);
        }
        else
        {
            flint_mpn_copyi(T, r, rn);
            tn = rn;
        }

        if (tn == 0)
        {
            /* v_c shares a factor with v_b */
            arg->bad = 1;
            break;
        }

        /* rc = (t * v_o) mod d */
        vo = pp + o * pl;
        von = C->level_len[k - 1][o];
        if (tn >= von)
            flint_mpn_mul(U, T, tn, vo, von);
        else
            flint_mpn_mul(U, vo, von, T, tn);
        un = tn + von;
        MPN_NORM(U, un);

        if (un > dn || (un == dn && mpn_cmp(U, d, dn) >= 0))
            mpn_tdiv_qr(Q, rc, 0, U, un, d, dn);
        else
            flint_mpn_copyi(rc, U, un);
    }

    flint_free(T);
}

typedef struct
{
    flint_mpn_crt_struct * C;
    nn_srcptr R;
    slong rl;
    int bad;
}
_crt_mult_arg;

static void
_crt_mult_worker(slong i, void * varg)
{
    _crt_mult_arg * arg = (_crt_mult_arg *) varg;
    flint_mpn_crt_struct * C = arg->C;
    slong cl = C->crt_chunk_limbs;
    slong j0 = C->chunk_start[i];
    slong j1 = C->chunk_start[i + 1];
    nn_srcptr Mi = C->level_prod[0] + i * cl;
    slong Min = C->level_len[0][i];
    nn_srcptr R = arg->R + i * arg->rl;
    slong Rn = arg->rl;
    nn_ptr Mt = FLINT_ARRAY_ALLOC(cl + 2, ulong);
    slong j, l;

    MPN_NORM(R, Rn);

    for (j = j0; j < j1; j++)
    {
        ulong lj = C->leaf[j].mod.n, t, inv, r, w, cof;
        slong Mtn = Min;

        /* Mt = M_i / l_j, cof = (Mt^-1 mod l_j) (R mod l_j)^-1 mod l_j;
           the second inverse is the chunk cofactor w_i reduced mod l_j,
           and its failure detects a common factor with another chunk */
        mpn_divexact_1(Mt, Mi, Min, lj);
        MPN_NORM(Mt, Mtn);
        if (Mtn == 0)
            Mtn = 1;
        t = mpn_mod_1(Mt, Mtn, lj);
        if (n_gcdinv(&inv, t, lj) != 1)
            goto bad;
        r = (Rn == 0) ? 0 : mpn_mod_1(R, Rn, lj);
        if (n_gcdinv(&w, r, lj) != 1)
            goto bad;
        cof = nmod_mul(inv, w, C->leaf[j].mod);

        if (C->crt_leaf_mult != NULL)
        {
            flint_mpn_copyi(C->crt_leaf_mult + j * cl, Mt, Mtn);
            flint_mpn_zero(C->crt_leaf_mult + j * cl + Mtn, cl - Mtn);
        }

        /* per-prime coefficients: u_j = sum_t r_t coeff_t mod l_j */
        if (C->leaf_start[j + 1] == C->leaf_start[j] + 1)
        {
            l = C->leaf_start[j];
            C->prime_crt_coeff[l] = cof;
            if (C->prime_crt_shoup != NULL)
                C->prime_crt_shoup[l] = n_mulmod_precomp_shoup(cof, lj);
        }
        else
        {
            for (l = C->leaf_start[j]; l < C->leaf_start[j + 1]; l++)
            {
                ulong p = C->primes[l], q, qinv;
                q = lj / p;
                if (n_gcdinv(&qinv, q % p, p) != 1)
                    goto bad;
                C->prime_crt_coeff[l] = nmod_mul(q * qinv, cof, C->leaf[j].mod);
            }
        }

        /* single-prime leaves: full multipliers V = Mt cof mod M_i, so that
           no modular multiplication is needed at all */
        if (C->crt_mult != NULL)
        {
            nn_ptr V = C->crt_mult + j * cl;
            slong Vn;
            ulong q[2];

            flint_mpn_zero(V, cl);
            Mt[Mtn] = mpn_mul_1(Mt, Mt, Mtn, cof);
            Vn = Mtn + 1;
            MPN_NORM(Mt, Vn);
            if (Vn > Min || (Vn == Min && mpn_cmp(Mt, Mi, Min) >= 0))
                mpn_tdiv_qr(q, V, 0, Mt, Vn, Mi, Min);
            else
                flint_mpn_copyi(V, Mt, Vn);
        }
    }

    flint_free(Mt);
    return;

bad:
    flint_free(Mt);
    arg->bad = 1;
}

void flint_mpn_crt_init2(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes, int flags)
{
    flint_mpn_crt_init_tuned(C, primes, num_primes, flags,
        FLINT_MPN_CRT_CHUNK_BITS, FLINT_MPN_CRT_MOD_BASE_BITS,
        FLINT_MPN_CRT_PREINV_CUTOFF);
}

void flint_mpn_crt_init(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes)
{
    flint_mpn_crt_init2(C, primes, num_primes, FLINT_MPN_CRT_MOD | FLINT_MPN_CRT_CRT);
}

void flint_mpn_crt_init_tuned(flint_mpn_crt_t C, nn_srcptr primes,
        slong num_primes, int flags, slong crt_chunk_bits, slong mod_base_bits,
        slong preinv_cutoff)
{
    int want_mod = (flags & FLINT_MPN_CRT_MOD) != 0;
    int want_crt = (flags & FLINT_MPN_CRT_CRT) != 0;

    slong i, j, k, l, n, nl, nc, per, cnt, top;
    ulong prod, t;
    slong total_bits;

    if (num_primes < 1)
        flint_throw(FLINT_ERROR, "flint_mpn_crt_init: num_primes must be positive\n");

    if (!want_mod && !want_crt)
        flint_throw(FLINT_ERROR, "flint_mpn_crt_init: no operations selected\n");

    n = num_primes;
    C->flags = flags;
    C->num_primes = n;
    C->primes = FLINT_ARRAY_ALLOC(n, ulong);
    for (i = 0; i < n; i++)
    {
        if (primes[i] == 0)
            flint_throw(FLINT_ERROR, "flint_mpn_crt_init: moduli must be nonzero\n");
        C->primes[i] = primes[i];
    }

    /* leaves: batch consecutive primes whose product fits in a limb */
    C->leaf_start = FLINT_ARRAY_ALLOC(n + 1, slong);
    C->leaf = FLINT_ARRAY_ALLOC(n, flint_mpn_crt_leaf_struct);
    C->prime_mod = FLINT_ARRAY_ALLOC(n, nmod_t);
    C->prime_barrett = FLINT_ARRAY_ALLOC(n, ulong);

    for (i = 0; i < n; i++)
    {
        nmod_init(C->prime_mod + i, primes[i]);
        C->prime_barrett[i] = (primes[i] == 1) ? UWORD_MAX : n_barrett_precomp(primes[i]);
    }

    nl = 0;
    i = 0;
    total_bits = 0;
    while (i < n)
    {
        prod = primes[i];
        j = i + 1;
        while (j < n && j - i < FLINT_MPN_CRT_MAX_LEAF_PRIMES &&
               !n_mul_checked(&t, prod, primes[j]))
        {
            prod = t;
            j++;
        }
        C->leaf_start[nl] = i;
        {
            flint_mpn_crt_leaf_struct * L = C->leaf + nl;
            ulong q[4], a[4] = {0, 0, 1, 0};
            nmod_init(&L->mod, prod);
            L->r64 = (-prod) % prod;
            L->r128 = nmod_mul(L->r64, L->r64, L->mod);
            mpn_divrem_1(q, 0, a, 3, prod);
            L->qlo = q[0];
            L->qhi = q[1];
            L->npre = (prod == 1) ? UWORD_MAX : n_barrett_precomp(prod);
        }
        total_bits += FLINT_BITS - C->leaf[nl].mod.norm;
        nl++;
        i = j;
    }
    C->leaf_start[nl] = n;
    C->num_leaves = nl;
    C->all_single = (nl == n);

    /* chunks: a power of two number of chunks (so that the tree is
       balanced, which matters for the scaled remainder tree), leaves
       distributed as evenly as possible */
    nc = (total_bits + crt_chunk_bits - 1) / crt_chunk_bits;
    nc = FLINT_MAX(nc, 1);
    nc = UWORD(1) << FLINT_BIT_COUNT(nc - 1);   /* next power of two */
    while (nc > nl)
        nc /= 2;
    C->num_chunks = nc;
    C->chunk_start = FLINT_ARRAY_ALLOC(nc + 1, slong);
    for (i = 0; i <= nc; i++)
        C->chunk_start[i] = (i * nl) / nc;
    per = 0;
    for (i = 0; i < nc; i++)
        per = FLINT_MAX(per, C->chunk_start[i + 1] - C->chunk_start[i]);

    /* tree shape */
    C->num_levels = 1;
    for (cnt = nc; cnt > 1; cnt = (cnt + 1) / 2)
        C->num_levels++;
    C->level_count = FLINT_ARRAY_ALLOC(C->num_levels, slong);
    C->level_limbs = FLINT_ARRAY_ALLOC(C->num_levels, slong);
    C->level_prod = FLINT_ARRAY_ALLOC(C->num_levels, nn_ptr);
    C->level_len = FLINT_ARRAY_ALLOC(C->num_levels, slong *);
    C->level_use_preinv = FLINT_ARRAY_ALLOC(C->num_levels, int);
    C->level_prod_norm = FLINT_ARRAY_ALLOC(C->num_levels, nn_ptr);
    C->level_inv = FLINT_ARRAY_ALLOC(C->num_levels, nn_ptr);
    C->level_norm = FLINT_ARRAY_ALLOC(C->num_levels, ulong *);

    for (k = 0, cnt = nc; k < C->num_levels; k++, cnt = (cnt + 1) / 2)
    {
        C->level_count[k] = cnt;
        C->level_len[k] = FLINT_ARRAY_ALLOC(cnt, slong);
        C->level_use_preinv[k] = 0;
        C->level_prod_norm[k] = NULL;
        C->level_inv[k] = NULL;
        C->level_norm[k] = NULL;
    }

    /* level 0: chunk products. A chunk of per leaves has at most per limbs. */
    {
        nn_ptr tmp = FLINT_ARRAY_ALLOC(nc * (per + 1), ulong);
        slong maxlen = 1;

        for (i = 0; i < nc; i++)
        {
            nn_ptr v = tmp + i * (per + 1);
            slong vn = 1;
            v[0] = 1;
            for (j = C->chunk_start[i]; j < C->chunk_start[i + 1]; j++)
            {
                v[vn] = mpn_mul_1(v, v, vn, C->leaf[j].mod.n);
                vn += (v[vn] != 0);
            }
            C->level_len[0][i] = vn;
            maxlen = FLINT_MAX(maxlen, vn);
        }

        C->level_limbs[0] = maxlen;
        C->level_prod[0] = FLINT_ARRAY_ALLOC(nc * maxlen, ulong);
        for (i = 0; i < nc; i++)
        {
            flint_mpn_copyi(C->level_prod[0] + i * maxlen, tmp + i * (per + 1), C->level_len[0][i]);
            flint_mpn_zero(C->level_prod[0] + i * maxlen + C->level_len[0][i], maxlen - C->level_len[0][i]);
        }
        flint_free(tmp);
    }

    /* higher levels */
    for (k = 1; k < C->num_levels; k++)
    {
        slong pl = C->level_limbs[k - 1];
        slong maxlen = 1;
        nn_ptr tmp = FLINT_ARRAY_ALLOC(C->level_count[k] * (2 * pl), ulong);
        _crt_prod_arg arg;

        arg.C = C;
        arg.k = k;
        arg.tmp = tmp;
        _crt_parallel_do(_crt_prod_worker, &arg, C->level_count[k], pl);

        for (i = 0; i < C->level_count[k]; i++)
            maxlen = FLINT_MAX(maxlen, C->level_len[k][i]);

        C->level_limbs[k] = maxlen;
        C->level_prod[k] = FLINT_ARRAY_ALLOC(C->level_count[k] * maxlen, ulong);
        for (i = 0; i < C->level_count[k]; i++)
        {
            flint_mpn_copyi(C->level_prod[k] + i * maxlen, tmp + i * (2 * pl), C->level_len[k][i]);
            flint_mpn_zero(C->level_prod[k] + i * maxlen + C->level_len[k][i], maxlen - C->level_len[k][i]);
        }
        flint_free(tmp);
    }

    top = C->num_levels - 1;
    C->prod = C->level_prod[top];
    C->prod_len = C->level_len[top][0];

    C->crt_mult = NULL;
    C->crt_leaf_mult = NULL;
    C->crt_fixed = 0;
    C->prime_crt_coeff = NULL;
    C->prime_crt_shoup = NULL;
    C->crt_use_shoup = 0;
    C->crt_fixed_m = 0;
    C->crt_small = NULL;
    C->crt_chunk_limbs = C->level_limbs[0];
    C->prod_half = NULL;

    if (want_crt)
    {
        C->prod_half = FLINT_ARRAY_ALLOC(C->prod_len, ulong);
        mpn_rshift(C->prod_half, C->prod, C->prod_len, 1);

        /*
            Cofactors. For a node with product v, let R = (P/v) mod v.
            The chunk cofactors are w_i = R_i^-1 mod M_i; the tree
            combination A = B v_c + C v_b then implicitly produces
            (P/v_a)^-1 x mod v_a at each node, and x itself at the top.
            R is computed for all chunks with a remainder tree:
            R_b = ((R_a mod v_b) v_c) mod v_b, starting from R = 1 at the
            top. This only requires single-limb-size inverses at the end.
        */
        {
            nn_ptr R, R2;
            slong rl, rl2;

            rl = 1;
            R = FLINT_ARRAY_ALLOC(1, ulong);
            R[0] = 1;

            for (k = top; k > 0; k--)
            {
                _crt_rem_arg arg;

                rl2 = C->level_limbs[k - 1];
                R2 = FLINT_ARRAY_ALLOC(C->level_count[k - 1] * rl2, ulong);

                arg.C = C;
                arg.k = k;
                arg.R = R;
                arg.rl = rl;
                arg.R2 = R2;
                arg.rl2 = rl2;
                arg.bad = 0;
                _crt_parallel_do(_crt_rem_worker, &arg, C->level_count[k], rl2);

                flint_free(R);
                R = R2;
                rl = rl2;

                if (arg.bad)
                {
                    flint_free(R);
                    goto bad_moduli;
                }
            }

            /* chunk level multipliers. The chunk cofactor w_i = R_i^-1 mod M_i
               is only needed modulo each leaf, where it equals
               (R_i mod l_j)^-1 mod l_j, so no big inverses are needed.
               Decide first whether the fixed-length path applies. */
            {
                _crt_mult_arg arg;
                int fixed = 0;

                if (nc == 1 && C->prod_len <= FLINT_MPN_CRT_FIXED_MAX)
                {
                    ulong hb = flint_clz(C->prod[C->prod_len - 1]);
                    ulong lb = FLINT_BIT_COUNT((ulong) nl);
                    ulong maxleaf = 0;

                    for (j = 0; j < nl; j++)
                        maxleaf = FLINT_MAX(maxleaf, C->leaf[j].mod.n);

                    /* the templates accumulate in independent 128-bit pairs, so we
                       need nl * 2^64 * max(u_j) < 2^128 as well as nl * P < 2^(64 N) */
                    fixed = (hb >= lb && !n_mul_checked(&t, (ulong) nl, maxleaf));
                }

                C->crt_fixed = fixed;
                C->prime_crt_coeff = FLINT_ARRAY_ALLOC(n, ulong);
                if (fixed || !C->all_single)
                    C->crt_leaf_mult = FLINT_ARRAY_ALLOC(nl * C->crt_chunk_limbs, ulong);
                if (!fixed && C->all_single)
                    C->crt_mult = FLINT_ARRAY_ALLOC(n * C->crt_chunk_limbs, ulong);
                if (fixed && C->all_single)
                {
                    C->crt_use_shoup = 1;
                    for (l = 0; l < n; l++)
                        if (C->prime_mod[l].norm == 0)
                            C->crt_use_shoup = 0;
                    if (C->crt_use_shoup)
                        C->prime_crt_shoup = FLINT_ARRAY_ALLOC(n, ulong);
                }

                arg.C = C;
                arg.R = R;
                arg.rl = (nc == 1) ? 1 : C->level_limbs[0];
                arg.bad = 0;
                _crt_parallel_do(_crt_mult_worker, &arg, nc, C->level_limbs[0]);
                if (arg.bad)
                {
                    flint_free(R);
                    goto bad_moduli;
                }

                if (fixed)
                {
                    C->crt_fixed_m = 1;
                    for (j = 0; j < nl; j++)
                    {
                        slong Mjn = C->prod_len;
                        nn_srcptr Mj = C->crt_leaf_mult + j * C->crt_chunk_limbs;
                        MPN_NORM(Mj, Mjn);
                        C->crt_fixed_m = FLINT_MAX(C->crt_fixed_m, Mjn);
                    }
                }
            }

            flint_free(R);
        }

        /* small-value shortcut structure (see multi_crt.c) */
        if (C->prod_len >= FLINT_MPN_CRT_SMALL_CHECK_LIMBS)
        {
            slong bits = 0, k = 0;

            /* first k primes with product of at most a few limbs (leaving
               room for the slack the fixed-length code needs) */
            while (k < n && bits + FLINT_BIT_COUNT(primes[k]) <= FLINT_BITS * FLINT_MPN_CRT_SMALL_CHECK_Q_LIMBS - 8)
                bits += FLINT_BIT_COUNT(primes[k++]);

            if (k >= 1)
            {
                C->crt_small = flint_malloc(sizeof(flint_mpn_crt_struct));
                flint_mpn_crt_init2(C->crt_small, primes, k, FLINT_MPN_CRT_CRT);

                /* without the fixed-length code (fullword moduli) the check
                   costs a division; only worth it for larger products */
                if (!C->crt_small->crt_fixed && C->prod_len < FLINT_MPN_CRT_SMALL_CHECK_LIMBS_GENERIC)
                {
                    flint_mpn_crt_clear(C->crt_small);
                    flint_free(C->crt_small);
                    C->crt_small = NULL;
                }
            }
        }
    }

    C->mod_pow = NULL;
    C->mod_base_level = top;
    C->mod_base_limbs = C->mod_pow_limbs = C->level_limbs[top];
    C->mod_pow_slack = 0;

    if (want_mod)
    {
        /* mod basecase level */
        C->mod_base_level = 0;
        while (C->mod_base_level < top &&
               C->level_limbs[C->mod_base_level] * FLINT_BITS < mod_base_bits)
            C->mod_base_level++;
        C->mod_base_limbs = C->level_limbs[C->mod_base_level];
        C->mod_pow_limbs = C->mod_base_limbs;
        if (C->mod_base_level == top)
            C->mod_pow_limbs += FLINT_MIN(FLINT_MAX(16, C->prod_len), 64);

        C->mod_pow = FLINT_ARRAY_ALLOC(nl * C->mod_pow_limbs, ulong);
        C->mod_pow_slack = 1;
        for (j = 0; j < nl; j++)
        {
            nmod_t mod = C->leaf[j].mod;
            ulong two64 = C->leaf[j].r64;
            nn_ptr pw = C->mod_pow + j * C->mod_pow_limbs;
            pw[0] = (mod.n == 1) ? 0 : 1;
            for (k = 1; k < C->mod_pow_limbs; k++)
                pw[k] = nmod_mul(pw[k - 1], two64, mod);
            if (n_mul_checked(&t, (ulong) C->mod_pow_limbs, mod.n))
                C->mod_pow_slack = 0;
        }

        /* division data: levels mod_base_level .. top-1 are targets of reductions */
        for (k = C->mod_base_level; k < top; k++)
        {
            slong cnt = C->level_count[k];
            slong ll = C->level_limbs[k];

            if (ll < preinv_cutoff)
                continue;

            C->level_use_preinv[k] = 1;
            C->level_prod_norm[k] = FLINT_ARRAY_ALLOC(cnt * ll, ulong);
            C->level_inv[k] = FLINT_ARRAY_ALLOC(cnt * ll, ulong);
            C->level_norm[k] = FLINT_ARRAY_ALLOC(cnt, ulong);

            for (i = 0; i < cnt; i++)
            {
                nn_srcptr d = C->level_prod[k] + i * ll;
                nn_ptr dn = C->level_prod_norm[k] + i * ll;
                slong len = C->level_len[k][i];
                ulong norm = flint_clz(d[len - 1]);

                C->level_norm[k][i] = norm;
                if (norm != 0)
                    mpn_lshift(dn, d, len, norm);
                else
                    flint_mpn_copyi(dn, d, len);
                flint_mpn_zero(dn + len, ll - len);
                flint_mpn_preinvn(C->level_inv[k] + i * ll, dn, len);
                flint_mpn_zero(C->level_inv[k] + i * ll + len, ll - len);
            }
        }

    }

    /* workspace: two ping-pong level buffers plus per-node pointer/length
       tables plus scratch for divisions and multiplications */
    {
        slong maxlev = 0, maxcnt = 0;
        for (k = 0; k < C->num_levels; k++)
        {
            maxlev = FLINT_MAX(maxlev, C->level_count[k] * (C->level_limbs[k] + FLINT_MPN_CRT_SLOT_EXTRA));
            maxcnt = FLINT_MAX(maxcnt, C->level_count[k]);
        }
        C->work_level_limbs = maxlev;
        C->work_max_count = maxcnt;
        C->tmp_limbs = 2 * maxlev + 4 * maxcnt + 6 * (C->prod_len + FLINT_MPN_CRT_SLOT_EXTRA) + 16;
    }

    return;

bad_moduli:
    flint_throw(FLINT_ERROR, "flint_mpn_crt_init: moduli are not pairwise coprime\n");
}

void flint_mpn_crt_clear(flint_mpn_crt_t C)
{
    slong k;

    flint_free(C->primes);
    flint_free(C->leaf_start);
    flint_free(C->leaf);
    flint_free(C->prime_mod);
    flint_free(C->prime_barrett);
    flint_free(C->crt_leaf_mult);
    if (C->crt_small != NULL)
    {
        flint_mpn_crt_clear(C->crt_small);
        flint_free(C->crt_small);
    }
    flint_free(C->prime_crt_coeff);
    flint_free(C->prime_crt_shoup);
    flint_free(C->chunk_start);
    flint_free(C->crt_mult);
    for (k = 0; k < C->num_levels; k++)
    {
        flint_free(C->level_prod[k]);
        flint_free(C->level_len[k]);
        flint_free(C->level_prod_norm[k]);
        flint_free(C->level_inv[k]);
        flint_free(C->level_norm[k]);
    }
    flint_free(C->level_count);
    flint_free(C->level_limbs);
    flint_free(C->level_prod);
    flint_free(C->level_len);
    flint_free(C->level_use_preinv);
    flint_free(C->level_prod_norm);
    flint_free(C->level_inv);
    flint_free(C->level_norm);
    flint_free(C->mod_pow);
    flint_free(C->prod_half);
}
