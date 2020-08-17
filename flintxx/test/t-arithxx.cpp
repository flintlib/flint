/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "arithxx.h"

#include "helpers.h"

using namespace flint;

void
test_stirling()
{
#define STIRLING(func, matfunc, n, k) \
{\
    tassert(func##_vec(n, k).size() == k); \
    fmpz_vecxx v1(func##_vec(n, k).evaluate() /* test temporary alloc */); \
    for(slong i = 0;i < v1.size();++i) \
        tassert(v1[i] == func(n, i)); \
    tassert(func##_vec_next(func##_vec(n, k), n+1).size() == k); \
    fmpz_vecxx v2(func##_vec_next(v1, n+1)); \
    for(slong i = 0;i < v2.size();++i) \
        tassert(v2[i] == func(n+1, i)); \
    fmpz_vecxx v3(func##_vec(n, n+1)); \
    fmpz_vecxx v4(func##_vec_next(v3, n+1)); \
    tassert(v4.size() == n+2 && v3.size() == n+1); \
    for(slong i = 0;i < v4.size();++i) \
        tassert(v4[i] == func(n+1, i)); \
    tassert(matfunc(n, k).rows() == n && matfunc(n, k).cols() == k); \
    fmpz_matxx M(matfunc(n, k).evaluate() /* test temporaries */); \
    for(slong i = 0;i < M.rows();++i) \
        for(slong j = 0;j < M.cols();++j) \
            tassert(M.at(i, j) == func(i, j)); \
}
STIRLING(stirling_number_1u, stirling_matrix_1u, 10, 6)
STIRLING(stirling_number_1, stirling_matrix_1, 10, 6)
STIRLING(stirling_number_2, stirling_matrix_2, 10, 6)
}

void
test_bell()
{
    tassert(bell_number_vec(10).size() == 10);
    fmpz_vecxx v(bell_number_vec(10));
    tassert(v == bell_number_vec_recursive(10));
    tassert(v == bell_number_vec_multi_mod(10));
    for(slong i = 0;i < v.size();++i)
        tassert(v[i] == bell_number(static_cast<unsigned>(i)));
    tassert(bell_number(10u) == bell_number_bsplit(10u));
    tassert(bell_number(10u) == bell_number_multi_mod(10u));

    nmodxx_ctx p(1031);
    tassert(bell_number_nmod(10u, p) == nmodxx::red(bell_number(10u), p));
    nmod_vecxx v2(bell_number_nmod_vec(10u, p));
    tassert(v2.size() == 10);
    for(slong i = 0;i < v2.size();++i)
        tassert(v2[i] == nmodxx::red(v[i], p));

    tassert(v2 == bell_number_nmod_vec_series(10u, p));
    tassert(v2 == bell_number_nmod_vec_recursive(10u, p));

    tassert(fmpzxx(2).pow(
                static_cast<unsigned> (bell_number_size(10u)))
            > bell_number(10u));
}

void
test_bernoulli()
{
    tassert(bernoulli_number(10u).den() == bernoulli_number_denom(10u));
    tassert(fmpqxx(2, 1u).pow(
                static_cast<slong> (bernoulli_number_size(10u)))
            > bernoulli_number(10u));
    fmpq_polyxx poly(bernoulli_polynomial(10u));
    tassert(poly.degree() == 10);
    for(slong i = 0;i < poly.length();++i)
        tassert(poly.get_coeff(i)
                == fmpqxx(bin(10u, (unsigned) i)*bernoulli_number(10u - (unsigned) i)));
    tassert(bernoulli_number_vec(10u).size() == 10);
    for(unsigned i = 0;i < 10;++i)
        tassert(bernoulli_number_vec(10u)[i] == bernoulli_number(i));
}

void
test_euler()
{
    tassert(fmpzxx(2).pow(
                static_cast<unsigned> (euler_number_size(10u)))
            > euler_number(10u));
    fmpq_polyxx poly(euler_polynomial(10u));
    tassert(poly.degree() == 10);
    tassert(poly(fmpqxx(1, 2u))*fmpqxx(2, 1u).pow(10)
            == fmpqxx(euler_number(10u).evaluate(), fmpzxx(1)));
    tassert(euler_number_vec(10u).size() == 10);
    for(unsigned i = 0;i < 10;++i)
        tassert(euler_number_vec(10u)[i] == euler_number(i));
}

void
test_legendre()
{
    const unsigned short N = 10;
    fmpq_polyxx f; f = "3  -1 0 1";
    f = f.pow(N);
    for(unsigned i = 0;i < N;++i)
        f = f.derivative();
    tassert(legendre_polynomial(N)
            == f * fmpqxx(fmpzxx(1), (fmpzxx(2).pow(N)*fac(N)).evaluate()));

    fmpq_polyxx x; x = "2  0 1";
    f = chebyshev_t_polynomial(N);
    f = f(cos_series(x, N));
    f.truncate(N);
    tassert(f == cos_series(N*x, N));
    tassert(chebyshev_u_polynomial(N)*(N + 1)
            == chebyshev_t_polynomial(N+1u).derivative());
}

void
test_multiplicative()
{
    fmpzxx p(1031);
    tassert(euler_phi(p) == p-1u);
    tassert(moebius_mu(p) == -1);
    tassert(divisor_sigma(p, 4u) == 1u + p.pow(4u));
    tassert(divisors(p).to_string() == "2  1 " + p.to_string());

    fmpz_polyxx q, one;
    q = "2  0 1";
    one = "1  1";
    fmpz_polyxx res(q * (one - q).pow(24u) * (one - q*q).pow(24u));
    for(unsigned i = 3;i < 10;++i)
        res *= (one - q.pow(i)).pow(24u);
    res.truncate(10);
    tassert(res == ramanujan_tau_series(10));
    for(int i = 0;i < 10;++i)
        tassert(ramanujan_tau_series(i+1).get_coeff(i)
                == ramanujan_tau(fmpzxx(i)));
}

void
test_polys()
{
    // just very basic tests ...
    ulong N = 1234;
    tassert(cyclotomic_polynomial(N).degree() == euler_phi(fmpzxx(N)));
    tassert(cos_minpoly(N).degree() == euler_phi(fmpzxx(N))/2);

    tassert(swinnerton_dyer_polynomial(8u).degree() == fmpzxx(2).pow(8u));
}

void
test_dedekind()
{
    frandxx state;
    fmpzxx h = fmpzxx::randtest_unsigned(state, 10);
    fmpzxx k = fmpzxx::randtest_unsigned(state, 10);
    tassert(dedekind_sum_naive(h, k) == dedekind_sum(h, k));
    k /= gcd(h, k);
    tassert(dedekind_sum_naive(h, k) == dedekind_sum(h, k));
}

void
test_number_of_partitions()
{
    unsigned short N = 15;
    nmodxx_ctx p(1031);
    fmpz_vecxx v1(number_of_partitions_vec(N));
    nmod_vecxx v2(number_of_partitions_nmod_vec(N, p));
    tassert(v1.size() == v2.size() && v1.size() == N);
    for(unsigned i = 0;i < N;++i)
    {
        tassert(v1[i] == number_of_partitions(i));
        tassert(nmodxx::red(v1[i], p) == v2[i]);
    }
}

void
test_sum_of_squares()
{
    unsigned short N = 15;
    fmpz_vecxx v(sum_of_squares_vec(7u, N));
    tassert(v.size() == N);
    for(unsigned i = 0;i < N;++i)
        tassert(v[i] == sum_of_squares(7u, fmpzxx(i)));
}

int
main()
{
    std::cout << "arithxx....";

    tassert(primorial(4) == 2*3);
    tassert(harmonic_number(3)
            == fmpqxx(1, 1u) + fmpqxx(1, 2u) + fmpqxx(1, 3u));

    test_stirling();
    test_bell();
    test_bernoulli();
    test_euler();
    test_legendre();
    test_multiplicative();
    test_polys();
    test_dedekind();
    test_number_of_partitions();
    test_sum_of_squares();

    tassert(landau_function_vec(9)[8] == 15);

    std::cout << "PASS\n";
    return 0;
}
