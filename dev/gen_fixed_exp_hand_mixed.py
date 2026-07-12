#!/usr/bin/env python3
"""Emit mixed-precision hand-style exp series evaluators for arbitrary
(n, r), in the shape of src/fixed/exp_rs_opt_hand.inc.

exp(x) = 1 + x + x^2 V_2,  V_k = C_k + x V_{k+1},  C_k = 1/k!

An error in V_k reaches the result multiplied by x^k, so with x < 2^-r
and an output ulp of 2^-64n, holding the series error near 2^-(64n-6),

    bits(V_k) = 64 n - 6 - r k

and each level is evaluated at ceil(bits/64) limbs.  Levels of one and
two limbs stay in registers (n_mulhi, the 2x1 and sloppy 2x2 forms);
from three limbs up they live in arrays and use flint_mpn_mulhigh_n,
with the dedicated 3x2 form where a 2-limb operand embeds with a zero
low limb.  x^2 uses flint_mpn_sqrhigh.
"""
import sys
from math import factorial, lgamma

LOG2 = 0.6931471805599453


def minN(n, r):
    N = 1
    while r * N + lgamma(N + 1) / LOG2 < 64 * n:
        N += 1
    return N


def limbs_for(n, r, k):
    bits = 64 * n - 6 - r * k
    if bits <= 0:
        return 1
    return max(1, -(-bits // 64))


def const_limbs(k, L):
    """1/k! as L limbs, low word first"""
    v = (1 << (64 * L)) // factorial(k)
    return [(v >> (64 * i)) & ((1 << 64) - 1) for i in range(L)]


def emit(o, n, r, name):
    N = minN(n, r)
    # limb count per level, clamped to n and made non-increasing inward
    L = {}
    for k in range(2, N):
        L[k] = min(n, limbs_for(n, r, k))
    for k in range(N - 1, 2, -1):          # inner levels never wider
        L[k - 1] = max(L[k - 1], L[k])

    o("/* n = %d, r = %d: N = %d terms; level widths %s */"
      % (n, r, N, " ".join("V%d:%d" % (k, L[k]) for k in range(2, N))))
    o("static void")
    o("%s(nn_ptr res, nn_srcptr x)" % name)
    o("{")
    maxarr = max(L.values())
    o("    ulong v, w1, w0, h, l, cy;")
    if maxarr >= 3:
        o("    ulong A[%d], B[%d], T[%d], x2[%d];" % (n, n, n, n))
    o("    ulong X2 = x[%d], X1 = x[%d];" % (n - 1, n - 2))
    o("")
    o("    (void) v; (void) w1; (void) w0; (void) h; (void) l;")
    o("")

    # innermost level: V_{N-1} = C_{N-1}
    k = N - 1
    cur = L[k]
    c = const_limbs(k, cur)
    if cur == 1:
        o("    v = UWORD(0x%016x);" % c[0])
    elif cur == 2:
        o("    w1 = UWORD(0x%016x); w0 = UWORD(0x%016x);" % (c[1], c[0]))
    else:
        for i in range(cur):
            o("    A[%d] = UWORD(0x%016x);" % (i, c[i]))

    # Horner inward -> outward
    for k in range(N - 2, 1, -1):
        prev = cur
        cur = L[k]
        c = const_limbs(k, cur)
        o("")
        o("    /* V_%d = C_%d + x V_%d  (%d x %d)" % (k, k, k + 1, cur, prev)
          + " */")
        if cur == 1:
            o("    v = UWORD(0x%016x) + n_mulhi(X2, v);" % c[0])
        elif cur == 2 and prev == 1:
            o("    _fixed_mulhi_2x1(&h, &l, X2, X1, v);")
            o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (c[1], c[0]))
        elif cur == 2 and prev == 2:
            o("    _fixed_mulhi_2x2_sloppy(&h, &l, X2, X1, w1, w0);")
            o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (c[1], c[0]))
        elif cur == 3 and prev == 2:
            o("    _fixed_mulhi_3x2(&A[2], &A[1], &A[0], x[%d], x[%d], x[%d],"
              % (n - 1, n - 2, n - 3))
            o("        w1, w0);")
            for i in range(3):
                o("    B[%d] = UWORD(0x%016x);" % (i, c[i]))
            o("    mpn_add_n(A, A, B, 3);")
        else:
            # general: widen the previous value into cur limbs, mulhigh
            if prev <= 2:
                o("    B[0] = %s;" % ("UWORD(0)" if prev == 1 else "w0"))
                if prev == 1:
                    for i in range(1, cur - 1):
                        o("    B[%d] = UWORD(0);" % i)
                    o("    B[%d] = v;" % (cur - 1))
                else:
                    for i in range(1, cur - 2):
                        o("    B[%d] = UWORD(0);" % i)
                    o("    B[%d] = w0; B[%d] = w1;" % (cur - 2, cur - 1))
                o("    flint_mpn_mulhigh_n(T, x + %d, B, %d);" % (n - cur, cur))
            else:
                if prev < cur:
                    for i in range(cur - prev):
                        o("    B[%d] = UWORD(0);" % i)
                    for i in range(prev):
                        o("    B[%d] = A[%d];" % (cur - prev + i, i))
                    o("    flint_mpn_mulhigh_n(T, x + %d, B, %d);"
                      % (n - cur, cur))
                else:
                    o("    flint_mpn_mulhigh_n(T, x + %d, A, %d);"
                      % (n - cur, cur))
            for i in range(cur):
                o("    B[%d] = UWORD(0x%016x);" % (i, c[i]))
            o("    mpn_add_n(A, T, B, %d);" % cur)

    # tail: res = x + x^2 V_2
    o("")
    o("    /* x^2 V_2 */")
    if cur <= 2:
        # V_2 in registers -> widen to n limbs
        o("    B[0] = UWORD(0);")
        if cur == 1:
            for i in range(1, n - 1):
                o("    B[%d] = UWORD(0);" % i)
            o("    B[%d] = v;" % (n - 1))
        else:
            for i in range(1, n - 2):
                o("    B[%d] = UWORD(0);" % i)
            o("    B[%d] = w0; B[%d] = w1;" % (n - 2, n - 1))
        src = "B"
    elif cur < n:
        for i in range(n - cur):
            o("    B[%d] = UWORD(0);" % i)
        for i in range(cur):
            o("    B[%d] = A[%d];" % (n - cur + i, i))
        src = "B"
    else:
        src = "A"
    o("    flint_mpn_sqrhigh(x2, x, %d);" % n)
    o("    flint_mpn_mulhigh_n(T, x2, %s, %d);" % (src, n))
    o("    cy = mpn_add_n(res, x, T, %d);" % n)
    o("    res[%d] = 1 + cy;" % n)
    o("}")
    o("")
    return N


def main():
    out = []
    o = out.append
    ent = []
    for n in (4, 5):
        for r in range(8, 41):
            if minN(n, r) > 21 or r > 64 * n - 16:
                continue
            nm = "hand_%d_%d" % (n, r)
            N = emit(o, n, r, nm)
            ent.append((nm, n, r, N))
    o("typedef void (*hand_fn)(nn_ptr, nn_srcptr);")
    o("typedef struct { hand_fn f; int n; int r; int N; } hand_entry;")
    o("static hand_entry hand_all[] = {")
    for nm, n, r, N in ent:
        o("    { %s, %d, %d, %d }," % (nm, n, r, N))
    o("};")
    o("static int hand_count = %d;" % len(ent))
    print("\n".join(out))


if __name__ == "__main__":
    main()
