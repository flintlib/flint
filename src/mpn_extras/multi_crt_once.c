/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "mpn_extras/impl.h"

/*
    One-shot multi CRT with O(N) memory.

    The algorithm is the same as in flint_mpn_multi_crt (products up,
    remainder tree down for the cofactors R = (P/v) mod v, combination
    A = A_b v_c + A_c v_b up), but the tree is traversed depth first:
    at a node we compute the products of the two halves, the cofactors of
    the children, then solve the children one after the other (or in
    parallel), freeing everything not needed. Only a constant number of
    values of the size of the current node are live along any path, so
    the total memory is O(N) rather than O(N log N).

    Subproducts are computed twice (once by the parent to derive the
    child cofactors, once by the child itself), which costs one extra
    product tree, cheap compared to the divisions in the remainder tree.

    A node value A for primes [a, b) with product v satisfies
    A = x (P/v)^-1 (mod v) and A < 2^(64 + log b + depth) v, so it fits
    in (b - a) + 2 limbs.
*/

#define ONCE_LEAF_PRIMES 48
#define ONCE_PROD_BASECASE 16

typedef struct
{
    nn_srcptr res;
    nn_srcptr primes;
    int bad;
}
_once_ctx;

/* v = product of primes[a, b), b > a; v must have b - a limbs; returns length */
static slong
_once_prod(nn_ptr v, nn_srcptr primes, slong a, slong b)
{
    slong vn, i;

    if (b - a <= ONCE_PROD_BASECASE)
    {
        v[0] = primes[a];
        vn = 1;
        for (i = a + 1; i < b; i++)
        {
            v[vn] = mpn_mul_1(v, v, vn, primes[i]);
            vn += (v[vn] != 0);
        }
        return vn;
    }
    else
    {
        slong m = a + (b - a) / 2;
        slong ln, rn;
        nn_ptr l, r;

        /* left half computed in place in v, then moved to a temporary
           only as large as needed */
        ln = _once_prod(v, primes, a, m);
        l = FLINT_ARRAY_ALLOC(ln, ulong);
        flint_mpn_copyi(l, v, ln);
        r = FLINT_ARRAY_ALLOC(b - m, ulong);
        rn = _once_prod(r, primes, m, b);

        if (ln >= rn)
            flint_mpn_mul(v, l, ln, r, rn);
        else
            flint_mpn_mul(v, r, rn, l, ln);
        vn = ln + rn;
        vn -= (v[vn - 1] == 0);

        flint_free(l);
        flint_free(r);
        return vn;
    }
}

/* r = ((R mod d) * o) mod d, given R with Rn limbs, d with dn limbs
   (normalised) and o with on limbs; r must have dn limbs. Returns the
   normalised length, or -1 if the moduli are not coprime. */
static slong
_once_cofactor(nn_ptr r, nn_srcptr R, slong Rn, nn_srcptr d, slong dn, nn_srcptr o, slong on)
{
    nn_ptr U, Q;
    slong tn, un, rn;

    /* t = R mod d, computed into r */
    if (Rn > dn || (Rn == dn && mpn_cmp(R, d, dn) >= 0))
    {
        Q = FLINT_ARRAY_ALLOC(Rn - dn + 1, ulong);
        mpn_tdiv_qr(Q, r, 0, R, Rn, d, dn);
        flint_free(Q);
        tn = dn;
        MPN_NORM(r, tn);
    }
    else
    {
        flint_mpn_copyi(r, R, Rn);
        tn = Rn;
    }

    if (tn == 0)
        return -1;

    /* u = t * o */
    U = FLINT_ARRAY_ALLOC(tn + on, ulong);
    if (tn >= on)
        flint_mpn_mul(U, r, tn, o, on);
    else
        flint_mpn_mul(U, o, on, r, tn);
    un = tn + on;
    MPN_NORM(U, un);

    /* r = u mod d */
    if (un > dn || (un == dn && mpn_cmp(U, d, dn) >= 0))
    {
        Q = FLINT_ARRAY_ALLOC(un - dn + 1, ulong);
        mpn_tdiv_qr(Q, r, 0, U, un, d, dn);
        flint_free(Q);
        rn = dn;
    }
    else
    {
        flint_mpn_copyi(r, U, un);
        rn = un;
    }

    MPN_NORM(r, rn);
    flint_free(U);
    return rn;
}

/* leaf: A = sum_t r_t V_t with V_t = (v/p_t) ((v/p_t)^-1 (R mod p_t)^-1 mod p_t) mod v */
static slong
_once_leaf(nn_ptr A, nn_ptr v, slong * vn, slong a, slong b, nn_srcptr R, slong Rn, _once_ctx * ctx)
{
    slong cl, t, An, Mtn, Vn;
    nn_ptr Mt, V;
    ulong hi = 0, lo = 0, cy, q[2];

    *vn = cl = _once_prod(v, ctx->primes, a, b);

    Mt = FLINT_ARRAY_ALLOC(cl + 2, ulong);
    V = FLINT_ARRAY_ALLOC(cl, ulong);

    flint_mpn_zero(A, cl);

    for (t = a; t < b; t++)
    {
        ulong p = ctx->primes[t], c, inv, w, f;
        nmod_t mod;

        mpn_divexact_1(Mt, v, cl, p);
        Mtn = cl;
        MPN_NORM(Mt, Mtn);
        if (Mtn == 0)
            Mtn = 1;

        c = mpn_mod_1(Mt, Mtn, p);
        if (n_gcdinv(&inv, c, p) != 1)
            goto bad;
        c = (Rn == 0) ? 0 : mpn_mod_1(R, Rn, p);
        if (n_gcdinv(&w, c, p) != 1)
            goto bad;
        nmod_init(&mod, p);
        f = nmod_mul(inv, w, mod);

        Mt[Mtn] = mpn_mul_1(Mt, Mt, Mtn, f);
        Vn = Mtn + 1;
        MPN_NORM(Mt, Vn);
        flint_mpn_zero(V, cl);
        if (Vn > cl || (Vn == cl && mpn_cmp(Mt, v, cl) >= 0))
            mpn_tdiv_qr(q, V, 0, Mt, Vn, v, cl);
        else
            flint_mpn_copyi(V, Mt, Vn);

        FLINT_ASSERT(ctx->res[t] < p);
        cy = mpn_addmul_1(A, V, cl, ctx->res[t]);
        add_ssaaaa(hi, lo, hi, lo, UWORD(0), cy);
    }

    A[cl] = lo;
    A[cl + 1] = hi;
    An = cl + 2;
    MPN_NORM(A, An);

    flint_free(Mt);
    flint_free(V);
    return An;

bad:
    ctx->bad = 1;
    flint_free(Mt);
    flint_free(V);
    return 0;
}

typedef struct
{
    nn_ptr A;
    slong An;
    nn_ptr v;
    slong vn;
    slong a;
    slong b;
    nn_srcptr R;
    slong Rn;
    int free_R;         /* whether this node owns R and frees it once the
                           children's cofactors have been computed */
    _once_ctx * ctx;
}
_once_arg;

/* solves the node into arg->A / arg->v (v must have b - a limbs,
   A must have b - a + 2 limbs) */
static void
_once_solve(void * varg)
{
    _once_arg * arg = (_once_arg *) varg;
    _once_ctx * ctx = arg->ctx;
    slong a = arg->a, b = arg->b;
    slong m, n = b - a;
    nn_ptr vb, vc, Rb, Rc, Ac, T;
    slong vbn, vcn, Rbn, Rcn, Abn, Acn, Tn, An;
    _once_arg argb, argc;
    thread_pool_handle * handles = NULL;
    slong num_threads = 0;

    if (ctx->bad)
    {
        if (arg->free_R)
            flint_free((nn_ptr) arg->R);
        return;
    }

    if (n <= ONCE_LEAF_PRIMES)
    {
        arg->An = _once_leaf(arg->A, arg->v, &arg->vn, a, b, arg->R, arg->Rn, ctx);
        if (arg->free_R)
            flint_free((nn_ptr) arg->R);
        return;
    }

    m = a + n / 2;

    /* products of the halves; v_b is computed into the output slot */
    vb = arg->v;
    vc = FLINT_ARRAY_ALLOC(b - m, ulong);
    vbn = _once_prod(vb, ctx->primes, a, m);
    vcn = _once_prod(vc, ctx->primes, m, b);

    /* cofactors of the children */
    Rb = FLINT_ARRAY_ALLOC(vbn, ulong);
    Rc = FLINT_ARRAY_ALLOC(vcn, ulong);
    Rbn = _once_cofactor(Rb, arg->R, arg->Rn, vb, vbn, vc, vcn);
    Rcn = _once_cofactor(Rc, arg->R, arg->Rn, vc, vcn, vb, vbn);

    /* the parent's cofactor is no longer needed */
    if (arg->free_R)
        flint_free((nn_ptr) arg->R);

    if (Rbn < 0 || Rcn < 0)
    {
        ctx->bad = 1;
        flint_free(vc);
        flint_free(Rb);
        flint_free(Rc);
        return;
    }

    /* the children recompute their products into vb and vc */
    Ac = FLINT_ARRAY_ALLOC(b - m + 2, ulong);

    argb.A = arg->A; argb.v = vb; argb.a = a; argb.b = m; argb.R = Rb; argb.Rn = Rbn; argb.free_R = 1; argb.ctx = ctx;
    argc.A = Ac; argc.v = vc; argc.a = m; argc.b = b; argc.R = Rc; argc.Rn = Rcn; argc.free_R = 1; argc.ctx = ctx;

    if (n >= FLINT_MPN_CRT_PARALLEL_LIMBS && flint_get_num_threads() > 1)
        num_threads = flint_request_threads(&handles, 2);

    if (num_threads >= 1)
    {
        thread_pool_wake(global_thread_pool, handles[0], 0, _once_solve, &argb);
        _once_solve(&argc);
        thread_pool_wait(global_thread_pool, handles[0]);
    }
    else
    {
        _once_solve(&argb);
        _once_solve(&argc);
    }

    if (handles != NULL)
        flint_give_back_threads(handles, num_threads);

    if (ctx->bad)
    {
        flint_free(vc);
        flint_free(Ac);
        return;
    }

    Abn = argb.An;
    Acn = argc.An;
    vbn = argb.vn;
    vcn = argc.vn;

    /* A = A_b v_c + A_c v_b, v = v_b v_c */
    T = FLINT_ARRAY_ALLOC(n + 3, ulong);

    if (Abn == 0)
    {
        Tn = 0;
    }
    else
    {
        if (Abn >= vcn)
            flint_mpn_mul(T, arg->A, Abn, vc, vcn);
        else
            flint_mpn_mul(T, vc, vcn, arg->A, Abn);
        Tn = Abn + vcn;
        MPN_NORM(T, Tn);
    }

    if (Acn == 0)
    {
        An = 0;
    }
    else
    {
        if (Acn >= vbn)
            flint_mpn_mul(arg->A, Ac, Acn, vb, vbn);
        else
            flint_mpn_mul(arg->A, vb, vbn, Ac, Acn);
        An = Acn + vbn;
        MPN_NORM(arg->A, An);
    }

    if (An == 0)
    {
        flint_mpn_copyi(arg->A, T, Tn);
        An = Tn;
    }
    else if (Tn != 0)
    {
        ulong cy;
        if (An >= Tn)
            cy = mpn_add(arg->A, arg->A, An, T, Tn);
        else
        {
            cy = mpn_add(arg->A, T, Tn, arg->A, An);
            An = Tn;
        }
        if (cy)
            arg->A[An++] = cy;
    }

    FLINT_ASSERT(An <= n + 2);
    arg->An = An;

    /* v = v_b v_c (v_b is in the output slot) */
    if (vbn >= vcn)
        flint_mpn_mul(T, vb, vbn, vc, vcn);
    else
        flint_mpn_mul(T, vc, vcn, vb, vbn);
    arg->vn = vbn + vcn;
    arg->vn -= (T[arg->vn - 1] == 0);
    flint_mpn_copyi(arg->v, T, arg->vn);

    flint_free(T);
    flint_free(vc);
    flint_free(Ac);
}

int
flint_mpn_multi_crt_once(nn_ptr out, slong * outn, nn_ptr prod, nn_srcptr res,
                         nn_srcptr primes, slong num_primes, int sign)
{
    _once_ctx ctx;
    _once_arg arg;
    nn_ptr A, v;
    slong An, vn, i;
    ulong one = 1;
    int negative = 0;

    for (i = 0; i < num_primes; i++)
        if (primes[i] == 0)
            flint_throw(FLINT_ERROR, "flint_mpn_multi_crt_once: moduli must be nonzero\n");

    ctx.res = res;
    ctx.primes = primes;
    ctx.bad = 0;

    A = FLINT_ARRAY_ALLOC(num_primes + 3, ulong);
    v = FLINT_ARRAY_ALLOC(num_primes, ulong);

    /* top cofactor: (P/P) mod P = 1 */
    arg.A = A;
    arg.v = v;
    arg.a = 0;
    arg.b = num_primes;
    arg.R = &one;
    arg.Rn = 1;
    arg.free_R = 0;
    arg.ctx = &ctx;
    _once_solve(&arg);

    if (ctx.bad)
    {
        flint_free(A);
        flint_free(v);
        flint_throw(FLINT_ERROR, "flint_mpn_multi_crt_once: moduli are not pairwise coprime\n");
    }

    An = arg.An;
    vn = arg.vn;

    FLINT_ASSERT(vn <= num_primes);
    FLINT_ASSERT(An <= vn + 2);

    /* final reduction mod P */
    if (An < vn || (An == vn && mpn_cmp(A, v, vn) < 0))
    {
        flint_mpn_copyi(out, A, An);
        flint_mpn_zero(out + An, vn - An);
    }
    else
    {
        ulong q[4];
        FLINT_ASSERT(An - vn + 1 <= 4);
        mpn_tdiv_qr(q, out, 0, A, An, v, vn);
    }

    if (sign)
    {
        /* negative iff 2 x > P */
        nn_ptr half = A;    /* reuse */
        mpn_rshift(half, v, vn, 1);
        if (mpn_cmp(out, half, vn) > 0)
        {
            mpn_sub_n(out, v, out, vn);
            negative = 1;
        }
    }

    *outn = vn;
    if (prod != NULL)
        flint_mpn_copyi(prod, v, vn);

    flint_free(A);
    flint_free(v);
    return negative;
}

void
flint_mpn_multi_mod_once(nn_ptr out, nn_srcptr x, slong xn, nn_srcptr primes, slong num_primes)
{
    flint_mpn_crt_t C;

    flint_mpn_crt_init2(C, primes, num_primes, FLINT_MPN_CRT_MOD);
    flint_mpn_multi_mod(out, x, xn, C, NULL);
    flint_mpn_crt_clear(C);
}
