/* Hand-written exp series for n = 1 (r = 12) and n = 2 (r = 16),
   the source dev/tune_mp_real.py carries into exp_opt_1.c and
   exp_opt_2.c -- these two predate the series generator and beat
   its output (N reduced to 4 and 7, plain Horner, sloppy high
   products; see dev/gen_mp_real_exp_rs_hard.py). */

#define MP_REAL_EXP_E2 UWORD(0x8000000000000000)
#define MP_REAL_EXP_E3 UWORD(0x2aaaaaaaaaaaaaaa)
#define MP_REAL_EXP_E4 UWORD(0x0aaaaaaaaaaaaaaa)
#define MP_REAL_EXP_E5 UWORD(0x0222222222222222)
#define MP_REAL_EXP_E6 UWORD(0x005b05b05b05b05b)
#define MP_REAL_EXP_E7 UWORD(0x000d00d00d00d00d)
#define MP_REAL_EXP_D2_HI UWORD(0x8000000000000000)
#define MP_REAL_EXP_D2_LO UWORD(0x0000000000000000)
#define MP_REAL_EXP_D3_HI UWORD(0x2aaaaaaaaaaaaaaa)
#define MP_REAL_EXP_D3_LO UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_EXP_G8  UWORD(0x0001a01a01a01a01)
#define MP_REAL_EXP_G9  UWORD(0x00002e3bc74aad8e)
#define MP_REAL_EXP_G10 UWORD(0x0000049f93edde27)
#define MP_REAL_EXP_H4_HI UWORD(0x0aaaaaaaaaaaaaaa)
#define MP_REAL_EXP_H4_LO UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_EXP_H5_HI UWORD(0x0222222222222222)
#define MP_REAL_EXP_H5_LO UWORD(0x2222222222222222)
#define MP_REAL_EXP_H6_HI UWORD(0x005b05b05b05b05b)
#define MP_REAL_EXP_H6_LO UWORD(0x05b05b05b05b05b0)
#define MP_REAL_EXP_H7_HI UWORD(0x000d00d00d00d00d)
#define MP_REAL_EXP_H7_LO UWORD(0x00d00d00d00d00d0)

/* exp(x), x < 2^-12, one fraction limb.  N = 4 terms: the dropped
   x^5/5! is about 2^-67, below the 2^-64 output ulp.  Written as
   1 + x + x^2 (C2 + x C3 + x^2 C4), reusing x^2 so that the two small
   terms are formed independently rather than in a serial chain. */
void
_mp_real_exp_rs_opt_1(nn_ptr res, nn_srcptr x)
{
    ulong X = x[0], x2, p, h, l;

    x2 = n_mulhi(X, X);
    p = MP_REAL_EXP_E2 + n_mulhi(X, MP_REAL_EXP_E3)
                     + n_mulhi(x2, MP_REAL_EXP_E4);
    p = n_mulhi(x2, p);

    add_ssaaaa(h, l, UWORD(0), X, UWORD(0), p);
    res[0] = l;
    res[1] = 1 + h;
}

/* exp(x), x < 2^-16, two fraction limbs.  N = 7 terms: the dropped
   x^8/8! is about 2^-143, below the 2^-128 output ulp.

   R = C4 + x (C5 + x (C6 + x C7)) needs only 58 bits, so it and the
   two Horner steps below it run at a single limb, on x's top limb
   alone; Q = C3 + x R and P = C2 + x Q then need two. */
void
_mp_real_exp_rs_opt_2(nn_ptr res, nn_srcptr x)
{
    ulong X1 = x[1], X0 = x[0];
    ulong u, s, r, q1, q0, p1, p0, h, l, x2h, x2l;

    /* inner tail at one limb */
    u = MP_REAL_EXP_E6 + n_mulhi(X1, MP_REAL_EXP_E7);
    s = MP_REAL_EXP_E5 + n_mulhi(X1, u);
    r = MP_REAL_EXP_E4 + n_mulhi(X1, s);

    /* Q = C3 + x R (R is a single limb) */
    _mp_real_mulhi_2x1(&h, &l, X1, X0, r);
    add_ssaaaa(q1, q0, MP_REAL_EXP_D3_HI, MP_REAL_EXP_D3_LO, h, l);

    /* P = C2 + x Q */
    _mp_real_mulhi_2x2_sloppy(&h, &l, X1, X0, q1, q0);
    add_ssaaaa(p1, p0, MP_REAL_EXP_D2_HI, MP_REAL_EXP_D2_LO, h, l);

    /* x^2 P, then the fraction x + x^2 P */
    _mp_real_sqrhi_2x2(&x2h, &x2l, X1, X0);
    _mp_real_mulhi_2x2_sloppy(&h, &l, x2h, x2l, p1, p0);
    add_ssaaaa(h, l, h, l, X1, X0);

    res[0] = l;
    res[1] = h;
    res[2] = 1;
}
