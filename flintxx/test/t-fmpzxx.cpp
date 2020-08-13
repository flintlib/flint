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

#include "fmpzxx.h"
#include "flintxx/test/helpers.h"

#if !HAVE_FAST_COMPILER
#warning "Some tests are disabled because your compiler is slow."
#endif

using namespace flint;

void
test_printing()
{
    fmpzxx a(31);

    tassert(a.to_string() == "31");
    tassert(a.to_string(2) == "11111");

    std::ostringstream oss;
    oss << a << '\n' << std::oct << a << '\n'
        << std::hex << a << '\n' << std::dec << a;
    tassert(oss.str() == "31\n37\n1f\n31");

    fmpzxx b(-15);
    tassert((a + b).to_string() == "16");

    char astr[] = "15";
    tassert(fmpzxx(astr) == 15);
}

void
test_order()
{
    fmpzxx a(0);
    fmpzxx b(1);
    fmpzxx c(0);
    fmpzxx d(-1);

#define TO(xzero, zero, one, mone) \
    tassert(xzero == zero); \
    tassert(xzero != one); \
    tassert(!(xzero == one)); \
    tassert(!(xzero != zero)); \
    tassert(xzero < one); \
    tassert(xzero <= one); \
    tassert(xzero <= zero); \
    tassert(!(xzero < zero)); \
    tassert(!(xzero <= mone)); \
    tassert(xzero > mone); \
    tassert(xzero >= mone); \
    tassert(xzero >= zero); \
    tassert(!(xzero > zero)); \
    tassert(!(xzero >= one));

    TO(a, c, b, d);
    TO(a, 0, 1, -1);
    TO(a, 0u, 1u, -1);
    TO(a, (signed short)(0), (signed short)(1), (signed short)(-1));
    TO(a, (unsigned short)(0), (unsigned short)(1), -1);
    TO(a, WORD(0), WORD(1), WORD(-1));
    TO(0, c, b, d);
    TO(0u, c, b, d);
    TO(WORD(0), c, b, d);
    TO((short)0, c, b, d);
    TO((unsigned short)0, c, b, d);
}

void
test_conversion()
{
    fmpzxx a(4);
    tassert(a.to<slong>() == 4);
    tassert(a.to<ulong>() == 4);
    tassert(a.to<double>() == 4.0);
    // NB: to_string is tested in test_printing
}

void
test_initialisation_assignment()
{
    fmpzxx a(4), b(WORD(4)), c(4u), d("4");
    fmpzxx e(a);
    fmpzxx f, g, h, i;
    f = 4;
    g = WORD(4);
    h = 4u;
    i = fmpzxx("4");
    tassert(a == b && a == c&& a == d && a == e && a == f && a == g && a == h
            && a == i);

    // test deep copying of fmpzxx with more than one digit
    a = fmpzxx("100000000000000000000");
    b = a;
    tassert(a._fmpz()[0] != b._fmpz()[0]);
    fmpzxx j(a);
    tassert(a._fmpz()[0] != j._fmpz()[0]);
    tassert(a == b && a == j);

    // just here to test our assumptions on data format
    tassert(c._fmpz()[0] == d._fmpz()[0]);

    // some more exotic "assignments"
    a.set_ui_smod(15, 16);
    tassert(a == -1);
    a.set_uiui(27, 15);
    tassert(a % (fmpzxx(1) << FLINT_BITS) == 15);
    tassert((a >> FLINT_BITS) == 27);
    b.neg_uiui(27, 15);
    tassert(b == -a);
}

void
test_arithmetic()
{
#define TAC(seven, three) \
    tassert(seven + three == 10); \
    tassert(seven * three == 21);

#define TA(seven, three) \
    TAC(seven, three); \
    tassert(seven - three == 4); \
    tassert(seven / three == 2); \
    tassert(seven % three == 1)

    TA(fmpzxx(7), fmpzxx(3));
    TA(fmpzxx(7), 3u);
    TAC(UWORD(7), fmpzxx(3));

    // test signed builtins
    tassert(-7 + fmpzxx(3) == -fmpzxx(4));
    tassert(fmpzxx(3) + 7 == fmpzxx(10));
    tassert(7 + fmpzxx(3) == fmpzxx(10));
    tassert(fmpzxx(3) + (-7) == -fmpzxx(4));
    tassert(fmpzxx(3) - 7 == -4);
    tassert(-7 * fmpzxx(3) == -21);
    tassert(fmpzxx(7) * (WORD(-3)) == -21);
    tassert(fmpzxx(21) / -3 == -7);

    // test composite arithmetic
    fmpzxx a(3), b(7);
    tassert(3*(a + b) - (b + (a - 4u)) + ((-(a - b)) % (b / 2)) == 25);

    // test unary minus
    tassert(-a == -3);

    // test assignment arithmetic
#define TAA(op, res) \
    { \
        fmpzxx tmp1(10), tmp2(10), tmp3(10); \
        fmpzxx three(3); \
        tmp1 op three; \
        tmp2 op 3u; \
        tmp3 op three*1; \
        tassert(tmp1 == res); \
        tassert(tmp2 == res); \
        tassert(tmp3 == res); \
    }
    TAA(+=, 13);
    TAA(*=, 30);
    TAA(/=, 3);
    TAA(%=, 1);

    // shifting
    tassert((fmpzxx(1) << 10) == 1024);
    tassert((fmpzxx(1024) >> 9) == 2);

    // binary logic
    tassert((fmpzxx(1) | fmpzxx(2)) == 3);
    tassert((fmpzxx(3) & fmpzxx(5)) == 1);
    tassert((fmpzxx(17) ^ fmpzxx(23)) == (UWORD(17) ^ UWORD(23)));
    tassert(~fmpzxx(17) == ~WORD(17));
}

void
test_functions()
{
    fmpzxx a(2);
    fmpzxx b(16);
    fmpzxx_srcref c(a);

    tassert(sgn(a) == 1);
    tassert(size(b) == 1);
    tassert(val2(b) == 4);
    tassert(bits(b) == 5);
    tassert(sizeinbase(b, 2) == 5);
    tassert(b.tstbit(0) == false);

    tassert(pow(a, 4u) == 16);
    tassert(root(b, 4) == 2);
    tassert(root(b, (unsigned short)4) == 2);
    tassert(sqrt(b) == 4);
    tassert(sqrt(a) == 1);
    tassert(rfac(a, 3u) == 2*3*4);
    tassert(a + fac(4u) == 2 + 4*3*2);
    tassert(fib(4u) == 3);
    tassert(bin(4u, 2u) == 6);
    tassert(abs(a) == a);
    tassert(gcd(a, b) == 2);
    tassert(lcm(a, b) == 16);
    fmpzxx p(31);
    tassert((invmod(b % p, p) * b) % p == 1);
    tassert((negmod(b % p, p) + b) % p == 0);

    fmpzxx mb(-b);
    fmpzxx three(3);
    tassert(cdiv_q(mb, three) == -5);
    tassert(fdiv_r(mb, three) == 2);
    tassert(tdiv_q(mb, three) == -5);
    tassert(cdiv_q(mb, 3) == -5);
    tassert(tdiv_q(mb, 3) == -5);
    tassert(cdiv_q(mb, 3u) == -5);
    tassert(tdiv_q(mb, 3u) == -5);
    tassert(tdiv_q_2exp(mb, 2u) == -4);
    tassert(fdiv_r_2exp(mb, 2u) == 0);
    tassert(divexact(b, a) == 8 && divexact(b, 2) == 8 && divexact(b, 2u) == 8);

    // check a composite expression
    tassert(2u + sqrt(a + a) == 4);

    // check immediate functions
    tassert(divisible(b, a + a));
    tassert(divisible(b, a + a + 2u) == false);
    tassert(divisible(a + a, (unsigned short)2));
    tassert(divisible(-a + b, 3) == false);
    tassert(divisible(c, a) == true);

    fmpzxx f(15);
    tassert(clog(f, a) == 4 && clog(f, 2u) == 4);
    tassert(flog(f, a) == 3 && flog(f, 2u) == 3);
    tassert(2.7 < dlog(f) && dlog(f) < 2.8);

    tassert(mul2(a, 15u, 16u) == a*15*16);
    tassert(divexact2(b, 2u, 2u) == 4);
    tassert(powm(a, 4u, fmpzxx(5)) == 1);
    tassert(powm(a, fmpzxx(4), fmpzxx(5)) == 1);
    tassert(mul_tdiv_q_2exp(a, fmpzxx(-3), 2u) == -1);
    tassert(mul_tdiv_q_2exp(a, -3, 2u) == -1);

    fmpzxx q(1), rem(2);
    tassert(fdiv_qr(b, a) == ltuple(8, 0));
    tassert(tdiv_qr(b, a) == ltuple(8, 0));

    bool worked;
    ltupleref(worked, q) = sqrtmod(fmpzxx(4), fmpzxx(5));
    tassert(worked && (q == 2 || q == 3));
    tassert(sqrtrem(fmpzxx(5)) == ltuple(2, 1));
    // TODO gcdinv, xgcd

    // check member functions
    slong exp = 0;
    double d = fmpzxx(3).get_d_2exp(exp);
    double r = d * (WORD(1) << exp);
    tassert(r >= 2.9 && r <= 3.1);
    fmpzxx one(1), mone(-1), zero(0), two(2);
    tassert(one.is_one() && one.is_pm1() && one.is_odd() &&
            !one.is_even() && !one.is_zero());
    tassert(!zero.is_one() && !zero.is_pm1() && zero.is_zero());
    tassert(mone.is_pm1());
    tassert(two.is_even() && !two.is_odd());
    tassert(one.is_square() && !two.is_square());
    tassert(jacobi(b, p) == 1);
    tassert(p.is_probabprime());
    tassert(p.is_prime_pseudosquare());

    tassert((3*b).remove(a) == ltuple(4, 3));

    a = 17;
    a.clrbit(0);
    tassert(a == 16);
    a.combit(1);
    tassert(a == 18);
    tassert(a.popcnt() == 2);


    frandxx rand;
    mp_bitcnt_t bits = 123; // some extra space in two words
    std::vector<mp_bitcnt_t> arr(bits / FLINT_BITS + 2);
    fmpzxx tostore = fmpzxx::randtest_unsigned(rand, bits);
    bit_pack(arr, bits, tostore);
    tassert(tostore == fmpzxx::bit_unpack_unsigned(arr, bits));
    tassert(tostore == fmpzxx::bit_unpack(arr, bits).get<1>());
}

void
test_member_functions()
{
    // NB: only a sample, since these are done by macros anyway
    fmpzxx a(2), b(3);
    tassert((a*b).divisible(b));
    tassert((a*b).divisible(3));
    tassert(b.flog(3u) == 1);
    tassert(a.clog(a) == 1);
    tassert(a.rfac(3u) == rfac(a, 3u));
    tassert(a.gcd(b) == gcd(a, b));
    tassert(a.lcm(b) == lcm(a, b));
    tassert(a.cdiv_q(b) == cdiv_q(a, b));
    tassert(a.mul2(3u, 4u) == mul2(a, 3u, 4u));
    tassert(a.sqrtmod(b) == sqrtmod(a, b));
    tassert(b.sqrt() == sqrt(b));
}

template<class T>
void assert_is_fmpzxx(const T&)
{
    tassert(traits::is_fmpzxx<T>::val);
}
struct newtype {typedef void data_ref_t; typedef void data_srcref_t;};
void
test_traits()
{
    tassert(traits::is_fmpzxx<fmpzxx>::val == true);
    tassert(traits::is_fmpzxx<int>::val == false);
    assert_is_fmpzxx(fmpzxx(1) + fmpzxx(2));
    // the following does not even compile:
    //tassert((traits::is_fmpzxx<fmpzxx_expression<
    //            operations::immediate, newtype> >::val == false));
}

template<class T>
unsigned count_temporaries2(const T&)
{
    return T::ev_traits_t::temp_rule_t::temporaries_t::len
         // this term is always zero, but causes compiler error
         // if we are not actually in the ternary case
         + T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER;
}

void
test_temporaries()
{
    fmpzxx a, b, c;
    tassert(count_temporaries(a + b) == 0);
    tassert(count_temporaries(a + b + c + a + b + c) == 1);
    tassert(count_temporaries(((a / c) + (b % a)) / ((b + c) + (c / a))) == 3);
    tassert(count_temporaries((a/b) + (a/c) + (b/c) + (c/b)) == 2);
    tassert(count_temporaries2((a*b) + (a*c) + (b*c) + (c*b)) == 1);

    // test a bug in evaluate_2 (if addmul is used on the right, this can
    // be done with two temporaries, else need three)
    tassert(count_temporaries(((a+b)+(a+c)) + (((a+c) + a*c))) == 2);
}

void
test_ternary()
{
    fmpzxx b(2), c(3), d(4);

#define T0 fac(4u)
#define T1 (b + b + b)
#define T2 (T1 + T1)
#define T3 (T2 + T2)
#define T4 (T3 + T3)

    // The inner struct is a trickery to get gcc to free some of its data
    // structures. It reduces the resident set by 50%, and compile time by 75%.
#define TT3(m1, m2, m3, ntemps) \
    do{ struct inner { static void doit() { \
    fmpzxx b(2), c(3), d(4); \
    tassert(count_temporaries2(m1 + m2*m3) == ntemps); \
    tassert(b + (m1 + m2*m3) == 2 + m1.to<slong>() + m2.to<slong>()*m3.to<slong>()); \
    tassert(count_temporaries2(m1 + m3*m2) == ntemps); \
    tassert(b + (m1 + m3*m2) == 2 + m1.to<slong>() + m2.to<slong>()*m3.to<slong>()); \
    \
    tassert(count_temporaries2(m2*m3 + m1) == ntemps); \
    tassert(b + (m2*m3 + m1) == 2 + m1.to<slong>() + m2.to<slong>()*m3.to<slong>()); \
    tassert(count_temporaries2(m1 + m3*m2) == ntemps); \
    tassert(b + (m3*m2 + m1) == 2 + m1.to<slong>() + m2.to<slong>()*m3.to<slong>()); \
    \
    tassert(count_temporaries2(m1 - m2*m3) == ntemps); \
    tassert(b + (m1 - m2*m3) == 2 + m1.to<slong>() - m2.to<slong>()*m3.to<slong>()); \
    tassert(count_temporaries2(m1 - m3*m2) == ntemps); \
    tassert(b + (m1 - m3*m2) == 2 + m1.to<slong>() - m2.to<slong>()*m3.to<slong>()); \
    } }; inner::doit();} while(0)
#define TT(m1, m2, ntemps) TT3(m1, m2, d, ntemps)

    TT(T0, c, 1);
    TT(T1, c, 1);
    TT(T2, c, 2);
    TT(T3, c, 3);

    TT(T0, T0, 2);
    TT(T0, T1, 2);
    TT(T0, T2, 2);
    TT(T0, T3, 3);

#if HAVE_FAST_COMPILER
    TT(T1, T0, 2);
    TT(T1, T1, 2);
    TT(T1, T2, 2);
    TT(T1, T3, 3);

    TT(T2, T0, 2);
    TT(T2, T1, 2);
    TT(T2, T2, 3);
    TT(T2, T3, 3);

    TT(T3, T0, 3);
    TT(T3, T1, 3);
    TT(T3, T2, 3);
    TT(T3, T3, 4);

    // NB: TT3 is symmetric in m2 and m3
#define TT6(m1, m2, m3, ntemps) \
    TT3(m1, m2, m3, ntemps); \
    TT3(m2, m1, m3, ntemps); \
    TT3(m3, m1, m2, ntemps);

    TT6(fac(2u), fac(3u), fac(4u), 3);
    TT6(T1, T2, T3, 3);
    TT6(T1, T2, T4, 4);
    TT6(T1, (d+d+d) /* T1' */, T4, 4);
    TT6(T0, fac(2u), T2, 3);
    TT6(T0, T1, (d+d+d), 3);
    TT6(T1, T3, (T2 + T1) /* T3' */, 3);
#endif
}

void
test_ternary_assigments()
{
    fmpzxx a(2), b(3), c(4);
    tassert((a += b*c) == 14);
    tassert(a == 14);
    tassert((a -= b*c) == 2);
    tassert(a == 2);

    tassert((a += (b+b)*c) == 26);
    tassert((a -= (b+b)*c) == 2);
    tassert((a += c*(b+b)) == 26);
    tassert((a -= c*(b+b)) == 2);

    tassert((a += (b+b)*(c+c)) == 50);
    tassert((a -= (b+b)*(c+c)) == 2);

    // Make sure that the general rule kicks
    tassert((a += 3*c) == 14);
    tassert((a -= 3*c) == 2);
}

bool
try_implicit_conversion(fmpzxx_ref, const fmpzxx_ref&,
        fmpzxx_srcref, const fmpzxx_srcref&, fmpzxx_srcref, const fmpzxx_srcref&)
{
    return true;
}
void
test_references()
{
    fmpzxx a(2), b(3);
    fmpzxx_ref ar(a);
    fmpzxx_srcref acr(a), bcr(b);

    tassert(ar == 2 && bcr == 3 && ar == acr && ar ==  a);

    // test assignments
    fmpzxx d(ar);
    a = 4;
    tassert(ar == 4 && acr == 4 && d == 2);
    ar = 2;
    tassert(a == 2);

    // test some arithmetic
    tassert(a + bcr == 5);
    tassert(acr + bcr == 5);
    tassert(ar + 3u == 5);
    tassert(a < bcr && ar < b && acr < bcr);
    tassert(rfac(acr, 1u) == 2);

    ar = bin(4u, 2u); tassert(acr == 6);
    tassert(acr.to<slong>() == WORD(6));
    tassert(ar.to_string() == "6");
    ar += b;
    ar += bcr;
    tassert(a == 12);

    // test conversion of reference types
    tassert(try_implicit_conversion(a, a, a, a, ar, ar));
    tassert((!traits::_is_convertible<fmpzxx, fmpzxx_ref>::val));
    tassert((!traits::_is_convertible<fmpzxx, fmpzxx_srcref>::val));
    tassert((!traits::_is_convertible<fmpzxx_ref, fmpzxx_srcref>::val));

    // test creation from C types
    fmpzxx_ref ar2 = fmpzxx_ref::make(a._fmpz());
    fmpzxx_srcref acr2 = fmpzxx_srcref::make(acr._fmpz());
    a = 7;
    tassert(ar2 == 7 && acr2 == 7);

    // test swapping
    a = 2;
    b = 3;
    swap(b, ar);
    tassert(a == 3 && b == 2);
    // make sure ADL is preferred over the general version in std
    using namespace std;
    swap(b, b);
}

void
test_randomisation()
{
    frandxx rand;
    tassert(abs(fmpzxx::randbits(rand, 5)) <= 31);
    tassert(abs(fmpzxx::randbits(rand, 5)) >= 16);
    tassert(abs(fmpzxx::randtest(rand, 5)) <= 31);
    tassert(fmpzxx::randtest_unsigned(rand, 5) <= 31);
    tassert(fmpzxx::randtest_unsigned(rand, 5) >= 0);
    tassert(fmpzxx::randtest_not_zero(rand, 5) != 0);
    fmpzxx thirty(30);
    tassert(fmpzxx::randm(rand, thirty) >= 0);
    tassert(fmpzxx::randm(rand, thirty) <= 29);
    tassert(fmpzxx::randtest_mod(rand, thirty) >= 0);
    tassert(fmpzxx::randtest_mod(rand, thirty) <= 29);
    tassert(fmpzxx::randtest_mod_signed(rand, thirty) > -15);
    tassert(fmpzxx::randtest_mod_signed(rand, thirty) <= 15);
}

void
test_factoring()
{
    fmpz_factorxx f;f = factor(-2*3*5);
    tassert(f.size() == 3);
    tassert(f.sign() == -1);
    tassert(f.p(0) == 2 && f.exp(0) == 1);
    tassert(f.p(1) == 3 && f.exp(1) == 1);
    tassert(f.p(2) == 5 && f.exp(2) == 1);
    tassert(f == factor(fmpzxx(-2*3*5)));
    tassert(f == fmpz_factorxx(f));
    tassert(f == factor_trial_range(fmpzxx(-2*3*5), 0u, 2u).get<1>());
    tassert(factor_trial_range(fmpzxx(-2*3*5), 0u, 2u).get<0>());

    tassert(f.expand() == -2*3*5);
    tassert(f.expand() == f.expand_iterative());
    tassert(f.expand() == f.expand_multiexp());

    if(0)
        print(f); // make sure this compiles
}

void
test_crt()
{
    frandxx rand;
    fmpzxx x = fmpzxx::randtest_unsigned(rand, 25);
    std::vector<mp_limb_t> primes;
    primes.push_back(1031);
    primes.push_back(1033);
    primes.push_back(1039);

    fmpz_combxx comb(primes);
    std::vector<mp_limb_t> residues(primes.size());
    multi_mod(residues, x, comb);

    fmpzxx prod(1);
    fmpzxx res;
    for(unsigned i = 0;i < primes.size();++i)
    {
        tassert(residues[i] == x % primes[i]);
        res = res.CRT(prod, residues[i], primes[i], false);
        prod *= primes[i];
    }
    tassert(res == x);

    res = multi_CRT(residues, comb, false);
    tassert(res == x);
}

int
main()
{
    std::cout << "fmpzxx....";

    test_printing();
    test_order();
    test_conversion();
    test_initialisation_assignment();
    test_arithmetic();
    test_functions();
    test_member_functions();
    test_traits();
    test_temporaries();
    test_ternary();
    test_ternary_assigments();
    test_references();
    test_randomisation();
    test_factoring();
    test_crt();
    test_print_read(fmpzxx(-17));

    // TODO test that certain things *don't* compile?
    // TODO test enable_all_fmpzxx

    std::cout << "PASS" << std::endl;
}
