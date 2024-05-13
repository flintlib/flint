#
#   Copyright (C) 2023 Fredrik Johansson
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

"""
Generic version of mul_n:

void flint_mpn_mul_n_basecase(nn_ptr res, nn_srcptr u, nn_srcptr v, slong n)
{
    ulong b, a;
    slong i;

    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);

    for (i = 2; i < n; i++)
        NN_DOTREV_S3_A3_1X1(b, a, res[i], 0, b, a, u, v, i + 1);

    for (i = n; i < 2 * n - 2; i++)
        NN_DOTREV_S3_A3_1X1(b, a, res[i], 0, b, a, u + i - n + 1, v + i - n + 1, 2 * n - i - 1);

    NN_ADDMUL_S2_A2_1X1(res[2 * n - 1], res[2 * n - 2], b, a, u[n - 1], v[n - 1]);
}

The compiler may refuse to unroll the nested loops, so we generate the code.

Schema for a general n x m multiply (here n = 7, m = 4):

    v0   u0 u1 u2 u3 u4 u5 u6 .
    v1      u0 u1 u2 u3 u4 u5 u6 .
    v2         u0 u1 u2 u3 u4 u5 u6 .
    v3            u0 u1 u2 u3 u4 u5 u6 .
"""

def mul1v(n):
    print("ulong flint_mpn_mul_%i_1v(nn_ptr res, nn_srcptr u, ulong v0)" % n)
    print("{")
    print("    ulong a;")
    print("    NN_MUL_1X1(a, res[0], u[0], v0);")
    for i in range(1, n-1):
        print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (i, i))
    print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (n - 1, n - 1))
    print("    return a;")
    print("}")

def mulnm(n, m):
    if m == 1:
        print("ulong flint_mpn_mul_%i_1(nn_ptr res, nn_srcptr u, nn_srcptr v)" % n)
        print("{")
        print("    ulong a, v0 = v[0];")
        print("    NN_MUL_1X1(a, res[0], u[0], v0);")
        for i in range(1, n-1):
            print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (i, i))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], 0, a, u[%i], v0);" % (n, n - 1, n - 1))
        print("    return res[%i];" % n)
        print("}")
    elif m == 2:
        print("ulong flint_mpn_mul_%i_%i(nn_ptr res, nn_srcptr u, nn_srcptr v)" % (n, m))
        print("{")
        print("    ulong b, a;")
        print("    ulong w[2];")
        print("    w[0] = v[0];")
        print("    w[1] = v[1];")
        print("    NN_MUL_1X1(a, res[0], u[0], w[0]);")
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);")
        for i in range(2, m):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u, w, %i);" % (i, i + 1))
        for i in range(m, n):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, w, %i);" % (i, i - m + 1, m))
        for i in range(n, n+m-2):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, w + %i, %i);" % (i, i - m + 1, i - n + 1, n + m - i - 1))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], w[%i]);" % (n + m - 1, n + m - 2, n - 1, m - 1))
        print("    return res[%i];" % (n + m - 1))
        print("}")
    else:
        print("ulong flint_mpn_mul_%i_%i(nn_ptr res, nn_srcptr u, nn_srcptr v)" % (n, m))
        print("{")
        print("    ulong b, a;")
        print("    NN_MUL_1X1(a, res[0], u[0], v[0]);")
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);")
        for i in range(2, m):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u, v, %i);" % (i, i + 1))
        for i in range(m, n):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v, %i);" % (i, i - m + 1, m))
        for i in range(n, n+m-2):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v + %i, %i);" % (i, i - m + 1, i - n + 1, n + m - i - 1))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], v[%i]);" % (n + m - 1, n + m - 2, n - 1, m - 1))
        print("    return res[%i];" % (n + m - 1))
        print("}")

for n in range(2, 8):
    for m in range(1, n+1):
        if n >= 8 and m >= 5:
            print("void flint_mpn_mul_%i_%i(nn_ptr res, nn_srcptr u, nn_srcptr v)" % (n, m))
            print("{")
            print("    __gmpn_mul_basecase(res, u, %i, v, %i);" % (n, m))
            print("}")
        else:
            mulnm(n, m)
        print()
    print()

for n in range(8, 17):
    mulnm(n, 1)
    print()
    #mulnm(n, 2)
    #print()

def mulhigh_n(n):
    print("ulong flint_mpn_mulhigh_%i_generic(nn_ptr res, nn_srcptr u, nn_srcptr v)" % n)
    print("{")
    print("    ulong b, a, low;")
    print("    NN_DOTREV_S3_1X1_HIGH(b, a, u, v, %i);" % (n - 1))
    print("    NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, u, v, %i);" % n)
    for i in range(n - 2):
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v + %i, %i);" % (i, i + 1, i + 1, n - i - 1))
    print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], v[%i]);" % (n - 1, n - 2, n - 1, n - 1))
    print("    return low;")
    print("}")

for n in range(3, 21):
    mulhigh_n(n)
    print()
