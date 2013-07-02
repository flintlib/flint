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

#ifndef CXX_PROTOTYPE_H
#define CXX_PROTOTYPE_H

#include <cstdlib>

// TODO
// * non-explicit ctos
// * contexts

#include "cxx/expression.h"
#include "cxx/stdmath.h"

#include "fmpz.h"
// also gets us FLINT_MIN and FLINT_MAX

namespace flint {
// "concrete" expression classes

template<class Operation, class Data>
class mpz_expression
    : public expression<derived_wrapper<mpz_expression>, Operation, Data>
{
public:
    mpz_expression() {}
    template<class T>
    explicit mpz_expression(const T& t) : mpz_expression::expression(t) {}

    template<class T>
    mpz_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit mpz_expression(const Data& d) : mpz_expression::expression(d) {}

    template<class D, class O, class Da>
    friend class expression;
};
typedef mpz_expression<operations::immediate, fmpz_t> mpz;

///////////////////////////////////////////////////////////////////
// HELPERS
///////////////////////////////////////////////////////////////////
namespace traits {
template<class T, class Enable = void>
struct is_mpz
    : mp::equal_types<typename T::evaluated_t, mpz> { };
template<class T>
struct is_mpz<T, typename mp::disable_if<traits::is_expression<T> >::type>
    : false_ { };
} // traits
namespace mp {
template<class Out, class T1, class T2 = void>
struct enable_all_mpz
    : mp::enable_if<mp::and_<traits::is_mpz<T1>, traits::is_mpz<T2> >, Out> { };
template<class Out, class T>
struct enable_all_mpz<Out, T, void>
    : mp::enable_if<traits::is_mpz<T>, Out> { };
} // mp

///////////////////////////////////////////////////////////////////
// RULES
///////////////////////////////////////////////////////////////////
namespace rules {
template<>
struct empty_initialization<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_init(v._data());
    }
};

template<>
struct destruction<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_clear(v._data());
    }
};

FLINT_DEFINE_DOIT(initialization, mpz, mpz,
        fmpz_init_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_unsigned_integer<T>,
        fmpz_init_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_signed_integer<T>,
        fmpz_init(to._data());fmpz_set_si(to._data(), from))

// TODO should expression automatically deduce this from the assignment
// implementation (or the other way round)?
template<int n>
struct initialization<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_init(target._data());
        // TODO what about different bases
        fmpz_set_str(target._data(), source, 10);
    }
};

FLINT_DEFINE_DOIT(assignment, mpz, mpz, fmpz_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_unsigned_integer<T>,
        fmpz_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_signed_integer<T>,
        fmpz_set_si(to._data(), from))

template<int n>
struct assignment<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_set_str(target._data(), const_cast<char*>(source), 10);
    }
};

template<>
struct cmp<mpz, mpz>
{
    static int get(const mpz& l, const mpz& r)
    {
        return fmpz_cmp(l._data(), r._data());
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_si(v._data(), t);
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_ui(v._data(), t);
    }
};

template<>
struct to_string<mpz>
{
    static std::string get(const mpz& v, int base)
    {
        char* str = fmpz_get_str(0, base, v._data());
        std::string res(str);
        std::free(str);
        return res;
    }
};

FLINT_DEFINE_GET(conversion, slong, mpz, fmpz_get_si(from._data()))
FLINT_DEFINE_GET(conversion, ulong, mpz, fmpz_get_ui(from._data()))
FLINT_DEFINE_GET(conversion, double, mpz, fmpz_get_d(from._data()))

FLINT_DEFINE_BINARY_EXPR(plus, mpz,
        fmpz_add(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(plus, mpz, traits::is_unsigned_integer<T>,
        fmpz_add_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(times, mpz,
        fmpz_mul(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_unsigned_integer<T>,
        fmpz_mul_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_signed_integer<T>,
        fmpz_mul_si(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(minus, mpz,
        fmpz_sub(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(minus, mpz, traits::is_unsigned_integer<T>,
        fmpz_sub_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(divided_by, mpz,
        fmpz_fdiv_q(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_unsigned_integer<T>,
        fmpz_fdiv_q_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_signed_integer<T>,
        fmpz_fdiv_q_si(to._data(), e1._data(), e2))

// TODO this interpretation of mod is not the same as for builtin types!
FLINT_DEFINE_BINARY_EXPR(modulo, mpz,
        fmpz_mod(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(modulo, mpz, traits::is_unsigned_integer<T>,
        fmpz_mod_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_UNARY_EXPR(negate, mpz, fmpz_neg(to._data(), from._data()))

// standard math functions (c/f stdmath.h)
FLINT_DEFINE_BINARY_EXPR_COND(pow, mpz, traits::is_unsigned_integer<T>,
        fmpz_pow_ui(to._data(), e1._data(), e2))
FLINT_DEFINE_BINARY_EXPR_COND(root, mpz, traits::fits_into_slong<T>,
        fmpz_root(to._data(), e1._data(), e2))
FLINT_DEFINE_UNARY_EXPR(sqrt, mpz, fmpz_sqrt(to._data(), from._data()))

namespace rdetail {
// Evaluation using ternary operators is actually surprisingly hard.
// Consider e.g. a + b*c. The number of temporaries needed for this depends
// on whether or not b, c are immediates, and on the numbers of temporaries
// needed for each non-immediate expression.

// This struct deals with the difficulties in whether b or c might be
// immediate.
template<class Left, class rigth1_t, class right2_t,
    bool bimm, bool cimm>
struct ternary_hhelper;

// Case where both are immediate.
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper<Left, right1_t, right2_t, true, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    static const unsigned norig = evl::temporaries_t::len;
    static const unsigned ntemps = FLINT_MAX(norig, 1);
    typedef typename mp::make_homogeneous_tuple<mpz*, ntemps>::type
        temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, mpz* res,
            const mpz*& toright1, const mpz*& toright2)
    {
        evl::doit(left._data(), mp::htuples::extract<norig>(temps), res);
        toright1 = &right1;
        toright2 = &right2;
    }
};

// If c is immediate but b is not, there are still two subcases.
// Let t1 be the number of temporaries needed to evaluate a, and
// t2 the number for b. If t1 >= t2, then we need to evaluate a first.
// Otherwise b.
// In any case, the number of temporaries is at least two (for the two return
// values), and generically equal to the maximum of t1 and t2. If however
// t1 == t2, then we need an additional temporary.
template<class Left, class right1_t, class right2_t, bool t1_ge_t2>
struct ternary_hhelper_1imm;
// Case where t1 >= t2
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper_1imm<Left, right1_t, right2_t, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef typename right1_t::ev_traits_t::temp_rule_t evr;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evr::temporaries_t::len;
    // t1 >= t2

    template<class Temps>
    static void doit(const Left& left, const right1_t& right1,
            Temps temps, mpz* res,
            const mpz*& toright1)
    {
        evl::doit(left._data(), mp::htuples::extract<t1>(temps), res);
        typename Temps::tail_t nores = mp::htuples::removeres(temps, res);
        evr::doit(right1._data(),
                mp::htuples::extract<t2>(nores), nores.head);
        toright1 = nores.head;
    }
};
// Case where t1 < t2
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper_1imm<Left, right1_t, right2_t, false>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef typename right1_t::ev_traits_t::temp_rule_t evr;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evr::temporaries_t::len;
    // t1 < t2

    template<class Temps>
    static void doit(const Left& left, const right1_t& right1,
            Temps temps, mpz* res,
            const mpz*& toright1)
    {
        typedef typename Temps::tail_t tail_t;
        tail_t nores = mp::htuples::removeres(temps, res);
        evr::doit(right1._data(),
                mp::htuples::extract<t2>(temps), nores.head);
        toright1 = nores.head;
        evl::doit(left._data(),
                mp::htuples::extract<t1>(tail_t(res, nores.tail)), res);
    }
};
// Case where c is immediate.
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper<Left, right1_t, right2_t, false, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef mp::evaluation_helper<right1_t> evhr1;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evhr1::temporaries_t::len;
    static const unsigned ntemps = FLINT_MAX(2, FLINT_MAX(t1, t2) + (t1 == t2));

    typedef ternary_hhelper_1imm<Left, right1_t, right1_t, t1 >= t2> thh1;
    typedef typename mp::make_homogeneous_tuple<mpz*, ntemps>::type
        temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, mpz* res,
            const mpz*& toright1, const mpz*& toright2)
    {
        toright2 = &right2;
        thh1::doit(left, right1, temps, res, toright1);
    }
};

// Case where b is immediate.
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper<Left, right1_t, right2_t, true, false>
{
    typedef ternary_hhelper<Left, right2_t, right1_t, false, true> thh;
    typedef typename thh::temporaries_t temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, mpz* res,
            const mpz*& toright1, const mpz*& toright2)
    {
        thh::doit(left, right2, right1, temps, res, toright2, toright1);
    }
};

// Case where neither is immediate.
template<class Left, class right1_t, class right2_t>
struct ternary_hhelper<Left, right1_t, right2_t, false, false>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef typename right1_t::ev_traits_t::temp_rule_t evr1;
    typedef typename right2_t::ev_traits_t::temp_rule_t evr2;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evr1::temporaries_t::len;
    static const unsigned t3 = evr2::temporaries_t::len;

    // m1, m2, m3 is t1, t2, t3 reordered s.t. m1 >= m2 >= m3
    static const unsigned m1 = FLINT_MAX(t1, FLINT_MAX(t2, t3));
    static const unsigned m3 = FLINT_MIN(t1, FLINT_MIN(t2, t3));
    static const unsigned m2 = t1 + t2 + t3 - m1 - m3;

    // The following is obtained by case analysis
    static const unsigned ntemps =
        (t1 == t2 && t2 == t3) ? FLINT_MAX(3, t1+2) : // all equal
        ((m1 > m2 && m2 > m3)  ? FLINT_MAX(3, m1) :   // all distinct
         (m1 == m2 ? FLINT_MAX(m1+1, 3)               // first two equal
                   : FLINT_MAX(m1, FLINT_MAX(m2+2, 3))));   // second two equal
    typedef typename mp::make_homogeneous_tuple<mpz*, ntemps>::type
        temporaries_t;

    struct resaccess
    {
        mpz* res;
        resaccess(mpz* r) : res(r) {};

        template<class Eval, class Temps, class Data>
        typename Temps::tail_t doit(const Data& d, Temps temps)
        {
            Eval::doit(d,
                    mp::htuples::extract<Eval::temporaries_t::len>(temps),
                    res);
            return mp::htuples::extract<Temps::len-1>(temps);
        }
    };
    struct toaccess
    {
        const mpz*& right;
        toaccess(const mpz*& r) : right(r) {};

        template<class Eval, class Temps, class Data>
        typename Temps::tail_t doit(const Data& d, Temps temps)
        {
            Eval::doit(d,
                    mp::htuples::extract<Eval::temporaries_t::len>(temps),
                    temps.head);
            right = temps.head;
            return temps.tail;
        }
    };

    struct doit_really
    {
        template<class E1, class E2, class E3, class A1, class A2, class A3>
        static void doit(const E1& e1, const E2& e2,
                const E3& e3, temporaries_t temps,
                A1 a1, A2 a2, A3 a3)
        {
            typedef typename E1::ev_traits_t::temp_rule_t ev1;
            typedef typename E2::ev_traits_t::temp_rule_t ev2;
            typedef typename E3::ev_traits_t::temp_rule_t ev3;
            a3.template doit<ev3>(e3._data(), 
                    a2.template doit<ev2>(e2._data(),
                        a1.template doit<ev1>(e1._data(), temps)));
        }
    };
    struct dont_doit
    {
        template<class E1, class E2, class E3, class A1, class A2, class A3>
        static void doit(const E1& e1, const E2& e2,
                const E3& e3, temporaries_t temps,
                A1 a1, A2 a2, A3 a3)
        {
        }
    };

    template<class E1, class E2, class E3, class A1, class A2, class A3>
    static void doit_sort(const E1& e1, const E2& e2,
            const E3& e3, temporaries_t temps,
            A1 a1, A2 a2, A3 a3)
    {
        typedef typename E1::ev_traits_t::temp_rule_t ev1;
        typedef typename E2::ev_traits_t::temp_rule_t ev2;
        typedef typename E3::ev_traits_t::temp_rule_t ev3;
        static const unsigned u1 = ev1::temporaries_t::len;
        static const unsigned u2 = ev2::temporaries_t::len;
        static const unsigned u3 = ev3::temporaries_t::len;

        if(u1 < u2)
            return doit_sort(e2, e1, e3, temps, a2, a1, a3);
        if(u2 < u3)
            return doit_sort(e1, e3, e2, temps, a1, a3, a2);

        // If we reach this point, u1 >= u2 >= u3.
        // However, even if this is not the case, the following line (and
        // everything it instantiates) still has to compile.
        mp::if_v<(u1 >= u2 && u2 >= u3), doit_really, dont_doit>::type::doit(
                e1, e2, e3, temps, a1, a2, a3);
    }

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, mpz* res,
            const mpz*& toright1, const mpz*& toright2)
    {
        // We re-order the temporaries in such a way that res is at the
        // very end. When evaluating things in the correct order, it is then
        // always correct to take temporaries from the front, and drop them
        // from the front.
        temporaries_t temps_reordered = mp::concat_tuple<
            typename temporaries_t::tail_t,
            mp::make_tuple<mpz*>::type>::doit(
                    mp::htuples::removeres(temps, res),
                    mp::make_tuple<mpz*>::make(res));
        doit_sort(left, right1, right2, temps_reordered,
                resaccess(res), toaccess(toright1), toaccess(toright2));
    }
};

template<class Left, class Right1, class Right2,
    class disable_op = void, class Enable = void>
struct ternary_helper { };
template<class Left, class Right1, class Right2, class disable_op>
struct ternary_helper<Left, Right1, Right2, disable_op,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Left>,
        traits::is_expression<typename traits::basetype<Right1>::type>,
        traits::is_expression<typename traits::basetype<Right1>::type> > >::type>
{
    typedef typename traits::basetype<Right1>::type right1_t;
    typedef typename traits::basetype<Right2>::type right2_t;

    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef mp::evaluation_helper<right1_t> evhr1;
    typedef mp::evaluation_helper<right2_t> evhr2;
    typedef mp::enable_if<mp::and_<
        traits::is_homogeneous_tuple<typename evl::temporaries_t, mpz*>,
        traits::is_homogeneous_tuple<typename evhr1::temporaries_t, mpz*>,
        traits::is_homogeneous_tuple<typename evhr2::temporaries_t, mpz*>,
        traits::is_homogeneous_tuple<
            typename mp::make_tuple<
                typename evl::return_t,
                typename evhr1::type,
                typename evhr2::type>::type, mpz>,
        mp::not_<mp::equal_types<typename Left::operation_t, disable_op> >
      > > enable;

    typedef ternary_hhelper<Left, right1_t, right2_t,
            traits::is_immediate<right1_t>::val,
            traits::is_immediate<right2_t>::val> inner;
    typedef typename inner::temporaries_t temporaries_t;

    // evaluate left into res, rigth1 and right2 to arbitrary location,
    // set toright1, toright2 to these locations
    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, mpz* res,
            const mpz*& toright1, const mpz*& toright2)
    {
        inner::doit(left, right1, right2, temps, res, toright1, toright2);
    }
};
} // rdetail

// a +- b*c
template<class Op, class Left, class Right1, class Right2>
struct evaluation<Op,
    tuple<Left, tuple<
        mpz_expression<operations::times,
            tuple<Right1, tuple<Right2, empty_tuple> > >,
        // NB: there is no particular reason to have the enable_if here,
        //     many other similar places would do
        typename mp::enable_if<mp::or_<
                mp::equal_types<Op, operations::plus>,
                mp::equal_types<Op, operations::minus> >,
            empty_tuple>::type> >,
    true, 1,
    typename rdetail::ternary_helper<Left, Right1, Right2>::enable::type>
{
    static const unsigned TERNARY_OP_MARKER = 0;

    typedef mpz return_t;
    typedef rdetail::ternary_helper<Left, Right1, Right2> th;
    typedef typename th::temporaries_t temporaries_t;
    typedef tuple<Left, tuple<
        mpz_expression<operations::times,
            tuple<Right1, tuple<Right2, empty_tuple> > >,
        empty_tuple> > data_t;
    static const bool is_add = mp::equal_types<Op, operations::plus>::val;

    static void doit(const data_t& input, temporaries_t temps, return_t* res)
    {
        const mpz* left = 0;
        const mpz* right = 0;
        th::doit(input.first(), input.second()._data().first(),
                input.second()._data().second(), temps, res, right, left);
        if(is_add)
            fmpz_addmul(res->_data(), left->_data(), right->_data());
        else
            fmpz_submul(res->_data(), left->_data(), right->_data());
    }
};

// b*c + a
template<class Right, class Left1, class Left2>
struct evaluation<operations::plus,
    tuple<mpz_expression<operations::times,
            tuple<Left1, tuple<Left2, empty_tuple> > >,
        tuple<Right, empty_tuple> >,
    true, 1,
    typename rdetail::ternary_helper<
        Right, Left1, Left2, operations::times>::enable::type>
{
    static const unsigned TERNARY_OP_MARKER = 0;

    typedef mpz return_t;
    typedef rdetail::ternary_helper<Right, Left1, Left2> th;
    typedef typename th::temporaries_t temporaries_t;
    typedef tuple<mpz_expression<operations::times,
            tuple<Left1, tuple<Left2, empty_tuple> > >,
        tuple<Right, empty_tuple> > data_t;

    static void doit(const data_t& input, temporaries_t temps, return_t* res)
    {
        const mpz* left = 0;
        const mpz* right = 0;
        th::doit(input.second(), input.first()._data().first(),
                input.first()._data().second(), temps, res, right, left);
        fmpz_addmul(res->_data(), left->_data(), right->_data());
    }
};
} // rules

// Assignment-arithmetic using ternaries.
namespace detail {
// XXX this basically reproduces the binary expression rule logic
template<class Right1, class Right2, class Enable = void>
struct ternary_assign_helper { };

// both immediate
template<class Right1, class Right2>
struct ternary_assign_helper<const Right1&, const Right2&,
    typename mp::enable_if<mp::and_<
        traits::is_immediate_expr<Right1>,
        traits::is_immediate_expr<Right2> > >::type>
{
    typedef mp::enable_if<mp::and_<
        traits::is_mpz<Right1>, traits::is_mpz<Right2> >, mpz&> enable;

    const mpz& right1;
    const mpz& right2;

    ternary_assign_helper(const Right1& r1, const Right2& r2)
        : right1(r1), right2(r2) {};
    const mpz& getleft() {return right1;}
    const mpz& getright() {return right2;}
};

// left immediate
template<class Right1, class Right2>
struct ternary_assign_helper<const Right1&, Right2,
    typename mp::enable_if<mp::and_<
        traits::is_immediate_expr<Right1>,
        traits::is_lazy_expr<Right2> > >::type>
{
    typedef mp::and_<traits::is_mpz<Right1>, traits::is_mpz<Right2> > cond;
    typedef mp::enable_if<cond, mpz&> enable;
    
    const mpz& right1;
    mpz right2;
    ternary_assign_helper(const Right1& r1, const Right2& r2)
        : right1(r1)
    {
        Right2::ev_traits_t::evaluate_into_fresh(right2, r2);
    }
    const mpz& getleft() {return right1;}
    const mpz& getright() {return right2;}
};

// right immediate
template<class Right1, class Right2>
struct ternary_assign_helper<Right1, const Right2&,
    typename mp::enable_if<mp::and_<
        traits::is_immediate_expr<Right2>,
        traits::is_lazy_expr<Right1> > >::type>
{
    typedef ternary_assign_helper<const Right2&, Right1> tah_t;
    typedef mp::enable_if<typename tah_t::cond, mpz&> enable;
    tah_t tah;
    ternary_assign_helper(const Right1& r1, const Right2& r2)
        : tah(r2, r1) {}
    const mpz& getleft() {return tah.getright();}
    const mpz& getright() {return tah.getleft();}
};

template<class Right1, class Right2>
struct ternary_assign_helper<Right1, Right2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Right2>,
        traits::is_lazy_expr<Right1> > >::type>
{
    // TODO this is literally copied from the binary expression rule
    typedef mp::enable_if<mp::and_<
        traits::is_mpz<Right1>, traits::is_mpz<Right2> >, mpz&> enable;

    typedef typename Right1::ev_traits_t::temp_rule_t ev1_t;
    typedef typename Right2::ev_traits_t::temp_rule_t ev2_t;
    typedef typename ev1_t::temporaries_t temporaries1_t;
    typedef typename ev2_t::temporaries_t temporaries2_t;
    typedef mp::merge_tuple<mp::make_tuple<mpz*>::type, temporaries2_t> merger2;
    typedef mp::merge_tuple<tuple<mpz*, typename merger2::type>,
                 temporaries1_t> merger1;
    typedef typename merger1::type temporaries_t;
    typedef mp::back_tuple<temporaries_t> back_t;

    typename back_t::type backing;
    temporaries_t temps;
    mpz* ret1;
    mpz* ret2;

    ternary_assign_helper(const Right1& r1, const Right2& r2)
    {
        back_t::init(temps, backing, 0);
        temporaries1_t temps1 = merger1::get_second(temps);
        temporaries2_t temps2 =
            merger2::get_second(merger1::get_first(temps).tail);
        ret1 = merger1::get_first(temps).head;
        ret2 = merger2::get_first(merger1::get_first(temps).tail).head;
        ev1_t::doit(r1._data(), temps1, ret1);
        ev2_t::doit(r2._data(), temps2, ret2);
    }

    const mpz& getright() {return *ret2;}
    const mpz& getleft() {return *ret1;}
};
} // detail

// a += b*c
template<class Right1, class Right2>
inline typename detail::ternary_assign_helper<Right1, Right2>::enable::type
operator+=(mpz& left, const mpz_expression<operations::times,
        tuple<Right1, tuple<Right2, empty_tuple> > >& other)
{
    detail::ternary_assign_helper<Right1, Right2> tah(
            other._data().first(), other._data().second());
    fmpz_addmul(left._data(), tah.getleft()._data(),
            tah.getright()._data());
    return left;
}

// a -= b*c
template<class Right1, class Right2>
inline typename detail::ternary_assign_helper<Right1, Right2>::enable::type
operator-=(mpz& left, const mpz_expression<operations::times,
        tuple<Right1, tuple<Right2, empty_tuple> > >& other)
{
    detail::ternary_assign_helper<Right1, Right2> tah(
            other._data().first(), other._data().second());
    fmpz_submul(left._data(), tah.getleft()._data(),
            tah.getright()._data());
    return left;
}


///////////////////////////////////////////////////////////////////////////
// FUNCTIONS
///////////////////////////////////////////////////////////////////////////

// These functions evaluate immediately and do not yield mpzs

template<class T1, class T2>
inline typename mp::enable_all_mpz<bool, T1, T2>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible(t1.evaluate()._data(), t2.evaluate()._data());
}
template<class T1, class T2>
inline typename mp::enable_if<mp::and_<
    traits::is_mpz<T1>, traits::fits_into_slong<T2> >, bool>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible_si(t1.evaluate()._data(), t2);
}

// These functions are evaluated lazily

FLINT_DEFINE_UNOP(fac)
FLINT_DEFINE_BINOP(rfac)
FLINT_DEFINE_BINOP(bin)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND(rfac, mpz, traits::is_unsigned_integer<T>,
        fmpz_rfac_ui(to._data(), e1._data(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(fac, mpz, traits::is_unsigned_integer<T>,
        fmpz_fac_ui(to._data(), from))

template<class T1, class T2>
struct binary_expression<
    T1,
    typename mp::enable_if<
        mp::and_<
            traits::is_unsigned_integer<T1>,
            traits::is_unsigned_integer<T2> >,
        operations::bin>::type,
    T2>
{
    typedef mpz return_t;
    static void doit(mpz& to, const T1& t1, const T2& t2)
    {
        fmpz_bin_uiui(to._data(), t1, t2);
    }
};
} // rules

// TODO many more functions

} // flint

#endif
