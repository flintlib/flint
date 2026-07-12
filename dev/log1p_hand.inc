/* hand-written fixed_log1p_opt_1 / _2 bodies, SPECIALIZED to
   r = 16 (the runtime-r branches of the originals for
   r >= 64 were pruned as unreachable).  dev/tune_fixed.py
   substitutes @FN@ and @SERIES@ and copies each body into
   log1p_opt_<n>.c on --pin 16 --emit, delimiting them by the
   section comments that open each one. */

/* n = 1: everything in registers.  With r <= 48 there is a single
   reduction window; each step is a funnel shift of P's top limbs, a
   branch-free compare-select, and masked updates of the deficit, the
   product and the (single-limb) logarithm accumulator.  The final
   division is a 3-by-2 limb schoolbook step, and at t < 2^-16 the
   atanh series collapses to t + t^3/3 (the t^5 term is below the
   64-bit ulp already for r = 16). */
void
@FN@(nn_ptr res, nn_srcptr x)
{
    const int r = 16;
    ulong p0 = 0, d0 = x[0], a0 = 0;
    ulong t0, z, w, q, rr, h, l, sn1, sn0, n2, n1, s1, s0, cy;
    slong i, nc;
    int k;


    _fixed_exp_logs_ensure(1, r);
    nc = _fixed_exp_logs_n;

    for (i = 1; i <= r + 1; i++)
    {
        slong ii = FLINT_MIN(i, (slong) r);
        ulong v = MPN_RIGHT_SHIFT_LOW(UWORD(1), p0, (int) ii);
        ulong lt = _fixed_exp_logs[ii * nc + (nc - 1)];
        ulong m = -(ulong) (v <= d0);

        d0 -= v & m;
        p0 += v & m;      /* no carry: P + v <= X < 2 */
        a0 += lt & m;     /* no carry: acc < log 2 */
    }

    /* t = d0 2^64 / S, S = X + P = (s1 : s0) with s1 in {2, 3}:
       one normalized 3-by-2 division step with two corrections */
    add_ssaaaa(s1, s0, UWORD(1), x[0], UWORD(1), p0);
    if (d0 == 0)
    {
        res[0] = a0;
        return;
    }
    k = flint_clz(s1);
    sn1 = (s1 << k) | (s0 >> (FLINT_BITS - k));
    sn0 = s0 << k;
    n2 = d0 >> (FLINT_BITS - k);
    n1 = d0 << k;
    udiv_qrnnd(q, rr, n2, n1, sn1);
    umul_ppmm(h, l, q, sn0);
    while (h > rr || (h == rr && l > 0))
    {
        q--;
        sub_ddmmss(h, l, h, l, UWORD(0), sn0);
        cy = rr + sn1;
        if (cy < rr)
            break;      /* remainder overflowed: q is now correct */
        rr = cy;
    }
    t0 = q;

    /* atanh(t) = t + t^3/3 + O(t^5), t < 2^-16 */
    umul_ppmm(z, l, t0, t0);
    umul_ppmm(w, l, t0, z);
    res[0] = a0 + 2 * (t0 + w / 3);
}

/* n = 2: two-limb register version.  Window 0 (i < 64) works on
   two-limb funnel shifts of P and branch-free two-limb
   compare-selects; window 1 (i >= 64) reduces to a single
   significant limb.  The division falls back to mpn_tdiv_qr and the
   atanh series runs four terms (enough down to r = 16). */
void
@FN@(nn_ptr res, nn_srcptr x)
{
    const int r = 16;
    ulong p1 = 0, p0 = 0, d1 = x[1], d0 = x[0], a1 = 0, a0 = 0;
    ulong t2[2], w[2], S[3], nd[4], rem[4];
    slong i, nc;


    _fixed_exp_logs_ensure(2, r);
    nc = _fixed_exp_logs_n;

#define STEP2(ii, v1, v0) \
    do { \
        nn_srcptr Lq = _fixed_exp_logs + (ii) * nc + (nc - 2); \
        ulong bw, e1, e0, m; \
        sub_dddmmmsss(bw, e1, e0, UWORD(0), d1, d0, \
            UWORD(0), v1, v0); \
        m = ~bw;                          /* accept iff no borrow */ \
        d1 = (d1 & bw) | (e1 & m); \
        d0 = (d0 & bw) | (e0 & m); \
        add_ssaaaa(p1, p0, p1, p0, (v1) & m, (v0) & m); \
        add_ssaaaa(a1, a0, a1, a0, Lq[1] & m, Lq[0] & m); \
    } while (0)

    for (i = 1; i <= (slong) r; i++)      /* r < FLINT_BITS: one window */
    {
        ulong v1 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1, (int) i);
        ulong v0 = MPN_RIGHT_SHIFT_LOW(p1, p0, (int) i);
        STEP2(i, v1, v0);
    }

    /* one extra step at i = r absorbs the truncation creep */
    {
        ulong v1 = MPN_RIGHT_SHIFT_LOW(UWORD(1), p1, (int) r);
        ulong v0 = MPN_RIGHT_SHIFT_LOW(p1, p0, (int) r);
        STEP2(r, v1, v0);
    }

#undef STEP2

    /* t = D 2^128 / (X + P) via mpn_tdiv_qr; S = X + P in [2, 4) */
    add_sssaaaaaa(S[2], S[1], S[0], UWORD(1), x[1], x[0],
        UWORD(1), p1, p0);
    nd[0] = 0;
    nd[1] = 0;
    nd[2] = d0;
    nd[3] = d1;
    t2[0] = t2[1] = 0;
    if (d1 != 0 || d0 != 0)
        mpn_tdiv_qr(t2, rem, 0, nd, 4, S, 3);

    /* generated atanh series for t < 2^-16 (single denominator
       division by truncated-inverse multiplication) */
    @SERIES@(w, t2);

    /* res = acc + 2 atanh(t) */
    add_ssaaaa(res[1], res[0], a1, a0,
        (w[1] << 1) | (w[0] >> (FLINT_BITS - 1)), w[0] << 1);
}
