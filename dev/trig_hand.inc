/* hand-written n = 1 half-angle body; dev/tune_fixed.py
   substitutes @FN@, @R@, @SERIES@ and copies it into
   trig_opt_1.c on --emit */

/* n = 1 without a single mpn call past the shared reduction: the
   rotation and series were already in registers, and here the whole
   reconstruction follows in the same A/B/C form as the generic tail --
   two umul_ppmm chains for TT and DE, three n_mulhi for A, B, C, and
   ONE udiv_qrnnd per requested reciprocal against a divisor normalized
   to 64 bits by a shift whose exponent moves into the final doubling.
   The roundings cost a few ulp against budgets of 6r + 128 and
   8r + 256.  Measured at roughly a third of the generic path at this
   size: mpn dispatch and buffer traffic dominated the call. */
void
@FN@(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x)
{
    slong num, nc;
    const int r = @R@;
    int e;
    ulong t[1], wx[2], wy[2];
    ulong U, X0, X1, Y0, ph, pl, qh, ql, n1, TT, D1, DE;
    ulong A, B, C, S, W, R, rr, v;
    slong * used;
    TMP_INIT;

    _fixed_atans_ensure(1, r);
    nc = _fixed_atans_n;

    TMP_START;
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    t[0] = x[0] >> 1;
    num = _fixed_bitwise_reduce(t, 1, r, 1, _fixed_atans, nc, used);
    _fixed_tan_rotate_tab[1](wx, wy, used, num);

    @SERIES@(t, t);          /* u = tan(t') */
    U = t[0];

    X0 = wx[0]; X1 = wx[1]; Y0 = wy[0];

    /* TT = wy + wx u < 0.64 (one limb) and DE = wx - wy u in
       (0.72, 1.17) (fraction limb + 0/1 unit limb), halved together
       when DE reaches one so that DE is normalized */
    umul_ppmm(ph, pl, X0, U);
    add_ssaaaa(n1, TT, UWORD(0), ph, UWORD(0), X1 ? U : UWORD(0));
    TT += Y0;                           /* n1 + carry stay zero */
    umul_ppmm(qh, ql, Y0, U);
    sub_ddmmss(D1, DE, X1, X0, UWORD(0), qh);
    (void) pl; (void) ql; (void) n1;
    if (D1)
    {
        DE = (UWORD(1) << (FLINT_BITS - 1)) | (DE >> 1);
        TT >>= 1;
    }

    /* A = TT DE, B = DE^2, C = TT^2: sin = 2A/(B + C),
       cos = 1 - 2C/(B + C), tan = 2A/(B - C), exactly as in the
       generic tail */
    A = n_mulhi(TT, DE);
    B = n_mulhi(DE, DE);
    C = n_mulhi(TT, TT);

    if (ysin != NULL || ycos != NULL)
    {
        ulong cy;
        add_ssaaaa(cy, S, UWORD(0), B, UWORD(0), C);
        if (cy)
        {
            S = (UWORD(1) << (FLINT_BITS - 1)) | (S >> 1);
            e = 1;
        }
        else if (!(S >> (FLINT_BITS - 1)))
        {
            S <<= 1;
            e = -1;
        }
        else
            e = 0;
        udiv_qrnnd(R, rr, ~UWORD(0) >> 1, ~UWORD(0), S);

        if (ysin != NULL)
        {
            v = n_mulhi(A, R);          /* sin = 2^(2-e) A R < 0.842 */
            ysin[1] = v >> (FLINT_BITS - 2 + e);
            ysin[0] = v << (2 - e);
        }
        if (ycos != NULL)
        {
            v = n_mulhi(C, R) << (2 - e);   /* 2C/S < 0.46 */
            ycos[0] = -v;
            ycos[1] = (v == 0);         /* cos = 1 exactly when v = 0 */
        }
    }

    if (ytan != NULL)
    {
        W = B - C;                      /* DE^2 (1 - t^2) > 0.175 */
        e = flint_clz(W);
        W <<= e;
        udiv_qrnnd(R, rr, ~UWORD(0) >> 1, ~UWORD(0), W);
        v = n_mulhi(A, R);              /* tan = 2^(2+e) A R < 1.56 */
        ytan[1] = v >> (FLINT_BITS - 2 - e);
        ytan[0] = v << (2 + e);
    }

    TMP_END;
}
