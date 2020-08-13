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

#include "fmpq_matxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpq_matxx A(3, 4);
    tassert(A.rows() == 3 && A.cols() == 4);
    tassert(A.at(0, 0).is_zero());
    A.at(0, 0) = 1;

    fmpq_matxx B(A);
    tassert(B.rows() == 3 && B.cols() == 4);
    tassert(B.at(0, 0).is_one());
    B.at(0, 0) = 0;
    tassert(A.at(0, 0).is_one());

    tassert(fmpq_matxx::zero(3, 4).is_zero());
    fmpq_matxx eye = fmpq_matxx::one(4, 4);
    for(slong i = 0;i < eye.rows();++i)
        for(slong j = 0;j < eye.cols();++j)
            tassert(eye.at(i, j) == fmpqxx::integer(int(i == j)));
}

void
test_assignment()
{
    frandxx state;
    fmpz_matxx A = fmpz_matxx::randtest(3, 4, state, 10);
    fmpq_matxx Aq(A.rows(), A.cols());
    Aq = A;
    tassert(Aq == fmpq_matxx::integer_matrix(A));
}

void
test_conversion()
{
    frandxx state;
    fmpz_matxx A = fmpz_matxx::randtest(3, 4, state, 10);
    fmpq_matxx Aq = fmpq_matxx::integer_matrix(A);
    tassert(Aq.rows() == A.rows() && Aq.cols() == A.cols());
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            tassert(Aq.at(i, j) == fmpqxx::integer(A.at(i, j)));
    tassert(A == fmpz_matxx::from_integral_fraction(Aq));
    Aq.at(0, 0) = fmpqxx::frac(1, 2);
    assert_exception(fmpz_matxx::from_integral_fraction(Aq));

    tassert(Aq.numden_entrywise().get<0>().rows() == A.rows()
            && Aq.numden_entrywise().get<0>().cols() == A.cols()
            && Aq.numden_entrywise().get<1>().rows() == A.rows()
            && Aq.numden_entrywise().get<1>().cols() == A.cols());

    tassert(Aq.numden_matwise().get<0>().rows() == A.rows()
            && Aq.numden_matwise().get<0>().cols() == A.cols());

    tassert(Aq.numden_rowwise().get<0>().rows() == A.rows()
            && Aq.numden_rowwise().get<0>().cols() == A.cols()
            && Aq.numden_rowwise().get<1>().size() == A.rows());

    tassert(Aq.numden_colwise().get<0>().rows() == A.rows()
            && Aq.numden_colwise().get<0>().cols() == A.cols()
            && Aq.numden_colwise().get<1>().size() == A.cols());

    fmpz_matxx den(A.rows(), A.cols());
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            den.at(i, j) = 1u;
    den.at(0, 0) = 2u;

    A.at(0, 0) = 1;
    tassert(Aq.numden_entrywise() == ltupleref(A, den));
    tassert(Aq.num_rowwise().at(0, 0) == 1);
    tassert(Aq.num_colwise().at(0, 0) == 1);
    tassert(Aq.numden_matwise().get<1>() == 2);

    fmpz_vecxx rowdens(A.rows());
    rowdens[0] = 2;
    for(slong i = 1;i < A.rows();++i)
        rowdens[i] = 1;
    for(slong i = 1;i < A.cols();++i)
        A.at(0, i) *= 2;
    tassert(Aq.numden_rowwise() == ltupleref(A, rowdens));
    tassert(Aq.numden_colwise().get<1>()[0] == 2);
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
void
test_arithmetic()
{
    fmpq_matxx A(10, 10);
    fmpq_matxx v(10, 1);
    for(int i = 0;i < 10;++i)
        v.at(i, 0) = i;

    fmpzxx two(2);
    tassert(transpose(v).rows() == 1);
    tassert(v.transpose().cols() == 10);
    tassert((two*v).rows() == 10);
    tassert((v*two).rows() == 10);
    tassert((v*transpose(v)).rows() == 10 && (v*transpose(v)).cols() == 10);

    tassert(!has_explicit_temporaries(trace(transpose(v))));
    tassert(!has_explicit_temporaries(trace(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(A + v*transpose(v)));
    tassert(!has_explicit_temporaries(trace((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(trace(v*transpose(v) + v*transpose(v))));
    tassert(!has_explicit_temporaries(v*transpose(v) + v*transpose(v)));

    tassert(trace(transpose(v)).is_zero());
    tassert(trace(A + v*transpose(v)) == fmpqxx(285, 1u));
    tassert(trace(v*transpose(v) + A) == fmpqxx(285, 1u));
    tassert(trace(v*transpose(v) + v*transpose(v)) == fmpqxx(2*285, 1u));
    tassert(trace((A+A)*(two + two)).is_zero());

    for(int i = 0;i < 10; ++i)
        for(int j = 0; j < 10; ++j)
            A.at(i, j) = i*j;
    tassert(A == v*transpose(v));
    tassert(A != transpose(v)*v);
    A.at(0, 0) = 15;
    tassert(A != v*transpose(v));

    A.at(0, 0) = 0;
    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) *= two;
    tassert(A == v*transpose(v) + v*transpose(v));
    tassert(A - v*transpose(v) == v*transpose(v));
    tassert(((-A) + A).is_zero());
    tassert((A + A).at(0, 0) == A.at(0, 0) + A.at(0, 0));

    tassert((A + A) == A*two);
    tassert((two*A) / two == A);

    frandxx state;
    A.set_randtest(state, 15);
    fmpz_matxx B(A.rows(), A.cols());
    B.set_randtest(state, 15);
    fmpq_matxx C(A.rows(), A.cols());
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            C.at(i, j).num() = B.at(i, j);
    tassert(C*A == B*A && A*C == A*B);
    tassert(C.mul_direct(A) == C*A && C.mul_cleared(A) == C*A);
}

void
test_functions()
{
    fmpq_matxx A(2, 3), B(2, 2), empty(0, 15);
    B.at(0, 0) = 1;
    tassert(A.is_zero() && !A.is_empty() && !A.is_square());
    tassert(!B.is_zero() == B.is_square());
    tassert(empty.is_zero() && empty.is_empty());

    // transpose tested in arithmetic
    // mul tested in arithmetic
    // trace tested in arithmetic

    tassert(hilbert_matrix(4, 6).rows() == 4);
    tassert(hilbert_matrix(4, 6).cols() == 6);
    A.set_hilbert_matrix();
    fmpq_matxx H(hilbert_matrix(2, 3));
    tassert(A == H);
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.rows();++j)
            tassert(A.at(i, j).num() == 1 && A.at(i, j).den() == i + j + 1);

    tassert(A.is_integral() == false);

    frandxx rand;
    fmpz_matxx Bp(B.rows(), B.cols());
    Bp.set_randtest(rand, 10);
    for(slong i = 0;i < B.rows();++i)
        for(slong j = 0;j < B.rows();++j)
            B.at(i, j) = fmpqxx(Bp.at(i, j), fmpzxx(1));
    tassert(B.det().num() == Bp.det() && B.det().den() == 1);

    B.at(0, 0) = fmpqxx::randtest_not_zero(rand, 10);
    B.at(0, 1) = 0;
    B.at(1, 0) = fmpqxx::randtest(rand, 10);
    B.at(1, 1) = fmpqxx::randtest_not_zero(rand, 10);

    tassert(B*B.solve_fraction_free(A) == A);
    tassert(B*B.solve_dixon(A) == A);
    fmpq_matxx eye(B.rows(), B.cols());
    for(slong i = 0;i < B.rows();++i)
        eye.at(i, i) = 1;
    tassert(B*B.inv() == eye);

    assert_exception(fmpq_matxx(2, 2).solve_fraction_free(A).evaluate());
    assert_exception(fmpq_matxx(2, 2).solve_dixon(A).evaluate());

    // make sure this compiles
    if(0)
        print(B);
}

void
test_extras()
{
    fmpq_matxx A(10, 10), B(10, 10);
    frandxx rand;
    A.set_randtest(rand, 15);
    B.set_randtest(rand, 15);
    A.at(0, 0) = B.at(0, 0) + fmpqxx(1, 1u);

    fmpq_matxx_srcref Asr(A);
    fmpq_matxx_ref Br(B);

    tassert((A + A) + (B + B) == (Asr + Asr) + (Br + Br));

    Br = Asr;
    tassert(A == B);

    fmpq_matxx C(Asr);
    tassert(C == A);
    C.at(0, 0) += fmpqxx(2, 1u);
    tassert(C != A);
}

void
test_randomisation()
{
    frandxx rand;
    fmpq_matxx A(2, 2);
    A.set_randbits(rand, 5);
    tassert(height(A.at(0, 0)) <= 31);
    A.set_randtest(rand, 5);
    tassert(height(A.at(0, 0)) <= 31);

    frandxx rand2, rand3;
    A.set_randbits(rand2, 5);
    tassert(A == fmpq_matxx::randbits(2, 2, rand3, 5));
    A.set_randtest(rand2, 5);
    tassert(A == fmpq_matxx::randtest(2, 2, rand3, 5));
}

void
test_reduction_reconstruction()
{
    fmpq_matxx A(4, 7);
    frandxx state;
    A.set_randtest(state, 4);
    fmpzxx M(UWORD(123457891));
    fmpz_matxx Ar = fmpz_matxx::reduce(A, M);
    tassert(Ar.rows() == A.rows() && Ar.cols() == A.cols());
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            tassert(Ar.at(i, j) == A.at(i, j) % M);

    tassert(A == fmpq_matxx::reconstruct(Ar, M));
    // TODO test exception
}

void
test_row_reduction()
{
    frandxx state;
    fmpq_matxx A = fmpq_matxx::randtest(5, 5, state, 15);
    slong rank1, rank2;
    permxx p1(5), p2(5);
    fmpq_matxx res1(A.rows(), A.cols()), res2(A.rows(), A.cols());

    fmpq_matxx B(A);
    tassert(B.pivot(1, 1, &p1) == fmpq_mat_pivot(p1._data(), A._mat(), 1, 1));
    tassert(A == B);

    ltupleref(rank1, res1) = rref(A);
    rank2 = fmpq_mat_rref(res2._mat(), A._mat());
    tassert(rank1 == rank2 && res1 == res2);

    tassert(rref_classical(A) == rref(A) && rref_fraction_free(A) == rref(A));
}

void
test_unified_access()
{
    fmpq_matxx A(2, 2);
    const fmpq_matxx& Ar = A;
    const fmpq_matxx_ref Ar2(A);
    Ar2.at(0, 0) = fmpqxx::one();
    tassert(Ar.at(0, 0).is_one());
}

int
main()
{
    std::cout << "fmpq_matxx....";

    test_init();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_extras();
    test_randomisation();
    test_reduction_reconstruction();
    test_row_reduction();
    test_unified_access();

    std::cout << "PASS" << std::endl;
    return 0;
}
