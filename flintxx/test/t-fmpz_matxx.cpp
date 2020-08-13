/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpz_matxx.h"
#include "fmpz_vecxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpz_matxx A(3, 4);
    tassert(A.rows() == 3 && A.cols() == 4);
    tassert(A.at(0, 0) == 0);
    A.at(0, 0) = 1;

    fmpz_matxx B(A);
    tassert(B.rows() == 3 && B.cols() == 4);
    tassert(B.at(0, 0) == 1);
    B.at(0, 0) = 0;
    tassert(A.at(0, 0) == 1);

    tassert(fmpz_matxx::zero(3, 4).is_zero());
    fmpz_matxx eye = fmpz_matxx::one(4, 4);
    for(slong i = 0;i < eye.rows();++i)
        for(slong j = 0;j < eye.cols();++j)
            tassert(eye.at(i, j) == int(i == j));
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
template<class T, class Expr>
bool compare_temporaries(const Expr&)
{
    return mp::equal_types<T,
             typename Expr::ev_traits_t::rule_t::temporaries_t>::val;
}
void
test_arithmetic()
{
    fmpz_matxx A(10, 10);
    fmpz_matxx v(10, 1);
    for(unsigned i = 0;i < 10;++i)
        v.at(i, 0) = i;

    tassert(transpose(v).rows() == 1);
    tassert(v.transpose().cols() == 10);
    tassert((2*v).rows() == 10);
    tassert((v*2).rows() == 10);
    tassert((v*transpose(v)).rows() == 10 && (v*transpose(v)).cols() == 10);
    tassert(mul_classical(v, transpose(v)).rows() == 10);
    tassert(mul_multi_mod(v, transpose(v)).cols() == 10);

    tassert(!has_explicit_temporaries(trace(transpose(v))));
    tassert(!has_explicit_temporaries(trace(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(A + v*transpose(v)));
    tassert(!has_explicit_temporaries(trace((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(trace(v*transpose(v) + v*transpose(v))));
    tassert(!has_explicit_temporaries(v*transpose(v) + v*transpose(v)));
    tassert((compare_temporaries<tuple<fmpzxx*, empty_tuple> >(
                    ((A+A)*(fmpzxx(1)+fmpzxx(1))))));

    tassert(trace(transpose(v)) == 0);
    tassert(trace(A + v*transpose(v)) == 285);
    tassert(trace(v*transpose(v) + A) == 285);
    tassert(trace(v*transpose(v) + v*transpose(v)) == 2*285);
    tassert(trace((A+A)*(fmpzxx(1) + fmpzxx(1))) == 0);

    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) = i*j;
    tassert(A == v*transpose(v));
    tassert(A != transpose(v)*v);
    A.at(0, 0) = 15;
    tassert(A != v*transpose(v));

    A.at(0, 0) = 0;
    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) *= 2;
    tassert(A == v*transpose(v) + v*transpose(v));
    tassert(A - v*transpose(v) == v*transpose(v));
    tassert(((-A) + A).is_zero());
    tassert((A + A).at(0, 0) == A.at(0, 0) + A.at(0, 0));

    tassert((A + A) == 2*A && A*2 == A*2u && fmpzxx(2)*A == 2u*A);
    tassert((2*A).divexact(2) == A);
    tassert((2*A).divexact(2u) == A);
    tassert((2*A).divexact(fmpzxx(2)) == A);
}

void
test_functions()
{
    fmpz_matxx A(2, 3), B(2, 2), empty(0, 15);
    B.at(0, 0) = 1;
    tassert(A.is_zero() && !A.is_empty() && !A.is_square());
    tassert(!B.is_zero() == B.is_square());
    tassert(empty.is_zero() && empty.is_empty());

    // transpose tested in arithmetic
    // mul tested in arithmetic
    // trace tested in arithmetic

    frandxx rand;
    A.set_randtest(rand, 10);
    B.set_randtest(rand, 10);
    tassert(B*A == mul_classical(B, A));
    tassert(B*A == mul_multi_mod(B, A));

    tassert(sqr(B) == B*B);
    tassert(B.sqr().sqr() == pow(B, 4u));

    B.set_randrank(rand, 1, 10);
    tassert(!inv(B).get<0>());

    B.set_randdet(rand, fmpzxx(2*3*5));
    tassert(B.det() == 2*3*5);
    fmpz_matxx Binv(2, 2); bool worked; fmpzxx d;
    ltupleref(worked, Binv, d) = inv(B);
    tassert(worked && d.divisible(fmpzxx(2*3*5)));
    fmpz_matxx eye(2, 2);eye.at(0, 0) = 1;eye.at(1, 1) = 1;
    tassert((Binv * B).divexact(d) == eye);

    B.set_randdet(rand, fmpzxx(105));
    tassert(B.det() == B.det_bareiss());
    tassert(B.det() == B.det_cofactor());
    tassert(abs(B.det()) <= B.det_bound());
    tassert(B.det().divisible(B.det_divisor()));
    tassert(B.det() == B.det_modular(true));
    tassert(B.det() == B.det_modular_accelerated(true));
    tassert(B.det() == B.det_modular_given_divisor(fmpzxx(1), true));

    tassert(B.charpoly().get_coeff(0) == B.det());
    tassert(charpoly(B).get_coeff(1) == -B.trace());
    tassert(charpoly(B).lead() == 1);

    A.set_randrank(rand, 2, 10);
    tassert(rank(A) == 2);

    fmpz_matxx X(2, 3);
    ltupleref(worked, X, d) = solve(B, A);
    tassert(worked == true && (B*X).divexact(d) == A);
    ltupleref(worked, X, d) = B.solve_fflu(A);
    tassert(worked == true && (B*X).divexact(d) == A);
    ltupleref(worked, X, d) = B.solve_cramer(A);
    tassert(worked == true && (B*X).divexact(d) == A);
    tassert(solve(B, A).get<1>() == X);

    slong nullity;fmpz_matxx C(3, 3);
    tassert(nullspace(A).get<1>().rows() == 3);
    tassert(nullspace(A).get<1>().cols() == 3);
    ltupleref(nullity, C) = nullspace(A);
    tassert(nullity == 3 - rank(A));
    tassert(C.rank() == nullity);
    tassert((A*C).is_zero());

    // TODO test solve_dixon, solve_bound
}

void
test_extras()
{
    fmpz_matxx A(10, 10), B(10, 10);
    frandxx rand;
    A.set_randtest(rand, 15);
    B.set_randtest(rand, 15);
    A.at(0, 0) = B.at(0, 0) + 1u;

    fmpz_matxx_srcref Asr(A);
    fmpz_matxx_ref Br(B);

    tassert((A + A) + (B + B) == (Asr + Asr) + (Br + Br));

    Br = Asr;
    tassert(A == B);

    fmpz_matxx C(Asr);
    tassert(C == A);
    C.at(0, 0) += 2u;
    tassert(C != A);
}

void
test_randomisation()
{
    frandxx rand;
    fmpz_matxx A = fmpz_matxx::randbits(2, 2, rand, 5);
    tassert(abs(A.at(0, 0)) <= 31 && abs(A.at(0, 0)) >= 16);
    A.set_randtest(rand, 5);
    tassert(abs(A.at(0, 0)) <= 31);
    fmpz_matxx::randtest(2, 2, rand, 5);

    fmpz_matxx B(2, 3);
    B.set_randintrel(rand, 5);
    tassert(abs(B.at(0, 0)) <= 31);

    A.set_randsimdioph(rand, 5, 6);
    tassert(A.at(0, 0) == 64 && abs(A.at(0, 1)) <= 31);
    tassert(A.at(1, 0) == 0  && A.at(1, 1) == 32);

    // TODO set_randntrulike, set_randntrulike2, set_randajtai

    fmpz_vecxx v(2);v[0] = 5;v[1] = 7;
    A.set_randpermdiag(rand, v);
    tassert(A.at(0, 0) + A.at(0, 1) + A.at(1, 0) + A.at(1, 1) == 5 + 7);

    A.set_randrank(rand, 1, 5);
    tassert(abs(A.at(0, 0)) <= 31 && A.rank() == 1);
    tassert(rank(fmpz_matxx::randrank(5, 6, rand, 3, 10)) == 3);
    A.apply_randops(rand, 17);
    tassert(abs(A.at(0, 0)) <= 31 && A.rank() == 1);

    A.set_randdet(rand, fmpzxx(17));
    tassert(det(A) == 17);
    tassert(fmpz_matxx::randdet(5, 5, rand, fmpzxx(123)).det() == 123);
}

void
test_row_reduction()
{
    frandxx state;
    fmpz_matxx A = fmpz_matxx::randtest(5, 5, state, 15);
    slong rank1, rank2;
    fmpzxx den1, den2;
    fmpz_matxx res1(A.rows(), A.cols()), res2(A.rows(), A.cols());

    tassert(find_pivot_any(A, 2, 4, 1)
            == fmpz_mat_find_pivot_any(A._mat(), 2, 4, 1));
    tassert(A.fflu(0, false).get<1>().rows() == A.rows());
    permxx p1(5), p2(5);
    ltupleref(rank1, res1, den1) = fflu(A, &p1);
    rank2 = fmpz_mat_fflu(res2._mat(), den2._fmpz(), p2._data(),
            A._mat(), false);
    tassert(rank1 == rank2 && res1 == res2 && p1 == p2 && den1 == den2);
    tassert(rank1 == A.fflu(0, false).get<0>());

    ltupleref(rank1, res1, den1) = rref(A);
    rank2 = fmpz_mat_rref(res2._mat(), den2._fmpz(), A._mat());
    tassert(rank1 == rank2 && res1 == res2 && den1 == den2);

    fmpz_matxx B(A);
    fmpzxx n(1031);
    A.set_rref_mod(n, &p1);
    fmpz_mat_rref_mod(p2._data(), B._mat(), n._fmpz());
    tassert(A == B && p1 == p2);
}

void
test_printing()
{
    frandxx rand;
    fmpz_matxx A = fmpz_matxx::randtest(2, 2, rand, 5);
    test_print_read(A);
    A.set_one();
    tassert_fprint_pretty(A, "[[1 0]\n[0 1]\n]");
}

int
main()
{
    std::cout << "fmpz_matxx....";

    test_init();
    test_arithmetic();
    test_functions();
    test_extras();
    test_randomisation();
    test_row_reduction();
    test_printing();

    std::cout << "PASS" << std::endl;
    return 0;
}
