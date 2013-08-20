/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

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
}

void
test_assignment()
{
    // TODO
}

void
test_conversion()
{
    // TODO
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
    //tassert(mul_classical(v, transpose(v)).rows() == 10);
    //tassert(mul_multi_mod(v, transpose(v)).cols() == 10);

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

    std::cout << "PASS" << std::endl;
    return 0;
}
