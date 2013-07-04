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

// This file contains helpers for evaluating expression templates.

#ifndef CXX_EVALUATION_TOOLS_H
#define CXX_EVALUATION_TOOLS_H

#include "flint.h" // FLINT_MAX and FLINT_MIN

#include "cxx/expression_traits.h"
#include "cxx/mp.h"
#include "cxx/rules.h"
#include "cxx/tuple.h"

namespace flint {
namespace mp {
// Find the highest-priority implemented evaluation rule, if any.
// TODO move to tools?
template<class Op, class Data,
    bool result_is_temporary, unsigned min_prio = 0>
struct find_evaluation
{
private:
    typedef rules::evaluation<Op, Data, result_is_temporary, 2> r2;
    typedef rules::evaluation<Op, Data, result_is_temporary, 1> r1;
    typedef rules::evaluation<Op, Data, result_is_temporary, 0> r0;

    typedef traits::is_implemented<r2> i2;
    typedef traits::is_implemented<r1> i1;
    typedef traits::is_implemented<r0> i0;

public:
    typedef typename mp::select<rules::UNIMPLEMENTED, // TODO
        mp::and_v<i2, min_prio <= 2>, r2,
        mp::and_v<i1, min_prio <= 1>, r1,
        mp::and_v<i0, min_prio <= 0>, r0
      >::type type;
};
} // mp

namespace tools {
namespace tdetail {
template<class T, class U>
struct cmp_invert
{
    static int get(const T& t, const U& u)
    {
        return -rules::cmp<U, T>::get(u, t);
    }
};
}

// A version of the cmp rule which tries both argument orders
template<class T, class U, class Enable = void>
struct symmetric_cmp
    : mp::if_<traits::is_implemented<rules::cmp<T, U> >,
          rules::cmp<T, U>,
          typename mp::if_<traits::is_implemented<rules::cmp<U, T> >,
              tdetail::cmp_invert<T, U>,
              rules::UNIMPLEMENTED
            >::type
        >::type { };


// Finding a subexpression of precsribed type
namespace tdetail {
template<class T, class Data, class Enable = void>
struct find_subexpr_helper2;

template<class T, class Expr, class Enable = void>
struct find_subexpr_helper
{
    typedef find_subexpr_helper2<T, typename Expr::data_t> fsh;
    static const bool val = fsh::val;
    static const T& get(const Expr& e) {return fsh::get(e._data());}
};

template<class T>
struct find_subexpr_helper<T, T>
{
    static const bool val = true;
    static const T& get(const T& t) {return t;}
};

template<class T, class Expr>
struct find_subexpr_helper<T, Expr,
    typename mp::enable_if<mp::and_<traits::is_immediate<Expr>,
        mp::not_<mp::equal_types<T, Expr> > > >::type>
{
    static const bool val = false;
};

template<class T, class Data, class Enable>
struct find_subexpr_helper2
{
    typedef find_subexpr_helper2<T, typename Data::tail_t> fsh;
    static const bool val = fsh::val;
    static const T& get(const Data& d) {return fsh::get(d.tail);}
};

template<class T, class Head, class Tail>
struct find_subexpr_helper2<T, tuple<Head, Tail>,
    typename mp::enable_if<find_subexpr_helper<T,
        typename traits::basetype<Head>::type> >::type>
{
    static const bool val = true;
    typedef typename traits::basetype<Head>::type head_t;
    typedef find_subexpr_helper<T, head_t> fsh;
    static const T& get(const tuple<Head, Tail>& d) {return fsh::get(d.head);}
};

template<class T>
struct find_subexpr_helper2<T, empty_tuple>
{
    static const bool val = false;
};
} // tdetail

template<class T, class Expr>
inline const T& find_subexpr(const Expr& e)
{
    return tdetail::find_subexpr_helper<T, Expr>::get(e);
}

template<class Expr, class T>
struct is_super_sub_expr : tdetail::find_subexpr_helper<T, Expr> { };


// A helper to invoke htuples::fill with instantiate_temporaries
namespace tdetail {
template<class Expr>
struct fill_tmps_helper
{
    const Expr& expr;
    fill_tmps_helper(const Expr& e) : expr(e) {};

    template<class T>
    T create() const {return rules::instantiate_temporaries<Expr, T>::get(expr);}
};
} // tdetail
template <class Expr>
tdetail::fill_tmps_helper<Expr> temporaries_filler(const Expr& e)
{
    return tdetail::fill_tmps_helper<Expr>(e);
}


// A helper to "evaluate" a single term, independend of whether or not it is
// actually an expression template
template<class T, class Enable = void>
struct evaluation_helper
{
    typedef typename traits::basetype<T>::type type;
    typedef typename traits::forwarding<type>::type ftype;
    static ftype get(const type& t) {return t;}

    typedef empty_tuple temporaries_t;
};

template<class T>
struct evaluation_helper<T,
    typename mp::enable_if<traits::is_lazy_expr<T> >::type>
{
    typedef typename T::evaluated_t type;
    static type get(const T& t) {return t.evaluate();}

    typedef typename T::ev_traits_t::temp_rule_t::temporaries_t temporaries_t;
};


///////////////////////////////////////////////////////////////////////////
// Helper to evaluate two terms
///////////////////////////////////////////////////////////////////////////
// Both Expr1 and Expr2 are allowed to be const references to immediates,
// or immediates (e.g. int), or lazy expressions. I.e. anything that is
// natually stored in the data slots of an expression.

template<class Expr1, class Expr2, class Enable = void>
struct evaluate_2;

// Case where both are immediate.
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_immediate<typename traits::basetype<Expr1>::type>,
        traits::is_immediate<typename traits::basetype<Expr2>::type> > >::type>
{
public:
    typedef typename evaluation_helper<Expr1>::type return1_t;
    typedef typename evaluation_helper<Expr2>::type return2_t;
    typedef typename evaluation_helper<Expr1>::ftype ref1_t;
    typedef typename evaluation_helper<Expr2>::ftype ref2_t;
    typedef ref1_t arg1_t;
    typedef ref2_t arg2_t;

private:
    ref1_t ref1;
    ref2_t ref2;

public:
    typedef empty_tuple temporaries_t;
    evaluate_2(temporaries_t temps, ref1_t e1, ref2_t e2)
        : ref1(e1), ref2(e2) {}

    ref1_t get1() {return ref1;}
    ref2_t get2() {return ref2;}
};

// Case where right is immediate
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_immediate<typename traits::basetype<Expr2>::type> > >::type>
{
private:
    typedef typename Expr1::ev_traits_t::temp_rule_t rule_t;
    typedef evaluation_helper<Expr2> eh2_t;

public:
    typedef typename rule_t::return_t return1_t;
    typedef const return1_t& ref1_t;
    typedef typename eh2_t::ftype ref2_t;
    typedef typename eh2_t::type return2_t;
    typedef const Expr1& arg1_t;
    typedef ref2_t arg2_t;

private:
    typedef typename rule_t::temporaries_t temporaries1_t;
    typedef mp::merge_tuple<
        typename mp::make_tuple<return1_t*>::type,
        temporaries1_t> merger;
    ref2_t ref2;
    return1_t* ret1;

public:
    typedef typename merger::type temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, ref2_t e2)
        : ref2(e2)
    {
        ret1 = merger::get_first(temps).head;
        rule_t::doit(e1._data(), merger::get_second(temps), ret1);
    }

    ref1_t get1() {return *ret1;}
    ref2_t get2() {return ref2;}
};

// Case where left is immediate
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr2>,
        traits::is_immediate<typename traits::basetype<Expr1>::type> > >::type>
{
private:
    typedef evaluate_2<Expr2, Expr1> ev2_t;
    ev2_t ev2;

public:
    typedef typename ev2_t::return1_t return2_t;
    typedef typename ev2_t::return2_t return1_t;
    typedef typename ev2_t::ref1_t ref2_t;
    typedef typename ev2_t::ref2_t ref1_t;
    typedef typename ev2_t::arg1_t arg2_t;
    typedef typename ev2_t::arg2_t arg1_t;
    typedef typename ev2_t::temporaries_t temporaries_t;

    evaluate_2(temporaries_t temps, ref1_t e1, const Expr2& e2)
        : ev2(temps, e2, e1) {};
    ref1_t get1() {return ev2.get2();}
    ref2_t get2() {return ev2.get1();}
};

// Case where neither is immediate
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_lazy_expr<Expr2> > >::type>
{
private:
    typedef typename Expr1::ev_traits_t::temp_rule_t rule1_t;
    typedef typename Expr2::ev_traits_t::temp_rule_t rule2_t;

public:
    typedef typename rule1_t::return_t return1_t;
    typedef typename rule2_t::return_t return2_t;
    typedef const return1_t& ref1_t;
    typedef const return2_t& ref2_t;
    typedef const Expr1& arg1_t;
    typedef const Expr2& arg2_t;

private:
    typedef typename rule1_t::temporaries_t temporaries1_t;
    typedef typename rule2_t::temporaries_t temporaries2_t;

    template<class E1, class E2, class Enable>
    friend struct evaluate_2;

    // We can either evaluate the Expr1 first and then Expr2, or the other
    // way round. We would like to choose the most efficient strategy.
    // Since we have no access to other metrics, we compare the number of
    // temporaries required (see typedef of doit below).

    struct doit_1
    {
        typedef mp::merge_tuple<typename mp::make_tuple<return2_t*>::type,
                     temporaries2_t> merger2;
        typedef mp::merge_tuple<tuple<return1_t*, typename merger2::type>,
                     temporaries1_t> merger1;
        typedef typename merger1::type temporaries_t;

        static void init(temporaries_t temps, const Expr1& e1, const Expr2& e2,
                return1_t*& ret1, return2_t*& ret2)
        {
            temporaries1_t temps1 = merger1::get_second(temps);
            temporaries2_t temps2 =
                merger2::get_second(merger1::get_first(temps).tail);
            ret1 = merger1::get_first(temps).head;
            ret2 =
                merger2::get_first(merger1::get_first(temps).tail).head;
            rule1_t::doit(e1._data(), temps1, ret1);
            rule2_t::doit(e2._data(), temps2, ret2);
        }
    };

    struct doit_2
    {
        typedef typename evaluate_2<Expr2, Expr1>::doit_1 doit_other;
        typedef typename doit_other::temporaries_t temporaries_t;

        static void init(temporaries_t temps, const Expr1& e1, const Expr2& e2,
                return1_t*& ret1, return2_t*& ret2)
        {
            doit_other::init(temps, e2, e1, ret2, ret1);
        }
    };

    typedef typename mp::if_v<
        (doit_1::temporaries_t::len <= doit_2::temporaries_t::len),
        doit_1, doit_2>::type doit;

    return1_t* ret1;
    return2_t* ret2;

public:
    typedef typename doit::temporaries_t temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, const Expr2& e2)
    {
        doit::init(temps, e1, e2, ret1, ret2);
    }

    ref1_t get1() {return *ret1;}
    ref2_t get2() {return *ret2;}
};

///////////////////////////////////////////////////////////////////////////
// Helper to evaluate three homogeneous terms
///////////////////////////////////////////////////////////////////////////
//
// Evaluation using ternary operators is actually surprisingly hard.
// Consider e.g. a + b*c. The number of temporaries needed for this depends
// on whether or not b, c are immediates, and on the numbers of temporaries
// needed for each non-immediate expression.
namespace tdetail {
// This struct deals with the difficulties in whether b or c might be
// immediate.
template<class T, class Left, class rigth1_t, class right2_t,
    bool bimm, bool cimm>
struct ternary_hhelper;
// To be specialised below.
} // tdetail

// The following struct can be used to simplify writing evaluation rules which
// use ternary operations (addmul, submul).
//
// In the situation of a + b*c, the optimization can be applied if
// - the result goes to a temporary (i.e. we can write to it prematurely)
// - a is not an immediate
// - a, b, c are of the same type, and addmul is available for this type
// If so, one needs to evaluate a into the return location and b, c into
// temporaries; after that addmul can be applied.
//
// The ternary_helper facilitates both the checking if we are in the right
// situation and the intermediate evaluations. Instantiate it with
// "T" being your ground type (for which addmul is implemented), "Left" the type
// of a, "Right1" the type of b and "Right2" the type of c.
// Then the member enable::type can be used in SFINAE situations to
// conditionally enable a template only if we are in the addmul situation.
// The member type "temporaries_t" and static member function "doit" can be used
// to evaluate the intermediate terms.
//
// It may sometimes be useful to preclude a certain type of expression for a.
// (E.g. one needs rules for both a + b*c and b*c + a, but then which of these
// applies to b*c + a*d?) To do this, pass the operation you want to exclude in
// "disable_op".
//
// NOTE: in the current implementation, ternary_helper only works with
// *homogeneous* expressions. These are defined to be expressions evaluating to
// type T, which only need temporaries of type T.
// This condition is included in the checks done by the enable member type.
//
template<class T, class Left, class Right1, class Right2,
    class disable_op = void, class Enable = void>
struct ternary_helper { };
template<class T, class Left, class Right1, class Right2, class disable_op>
struct ternary_helper<T, Left, Right1, Right2, disable_op,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Left>,
        traits::is_expression<typename traits::basetype<Right1>::type>,
        traits::is_expression<typename traits::basetype<Right1>::type> > >::type>
{
    typedef typename traits::basetype<Right1>::type right1_t;
    typedef typename traits::basetype<Right2>::type right2_t;

    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef tools::evaluation_helper<right1_t> evhr1;
    typedef tools::evaluation_helper<right2_t> evhr2;
    typedef mp::enable_if<mp::and_<
        traits::is_homogeneous_tuple<typename evl::temporaries_t, T*>,
        traits::is_homogeneous_tuple<typename evhr1::temporaries_t, T*>,
        traits::is_homogeneous_tuple<typename evhr2::temporaries_t, T*>,
        traits::is_homogeneous_tuple<
            typename mp::make_tuple<
                typename evl::return_t,
                typename evhr1::type,
                typename evhr2::type>::type, T>,
        mp::not_<mp::equal_types<typename Left::operation_t, disable_op> >
      > > enable;

    typedef tdetail::ternary_hhelper<T, Left, right1_t, right2_t,
            traits::is_immediate<right1_t>::val,
            traits::is_immediate<right2_t>::val> inner;
    typedef typename inner::temporaries_t temporaries_t;

    // evaluate left into res, rigth1 and right2 to arbitrary location,
    // set toright1, toright2 to these locations
    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, T* res,
            const T*& toright1, const T*& toright2)
    {
        inner::doit(left, right1, right2, temps, res, toright1, toright2);
    }
};

namespace tdetail {
// Case where both are immediate.
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper<T, Left, right1_t, right2_t, true, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    static const unsigned norig = evl::temporaries_t::len;
    static const unsigned ntemps = FLINT_MAX(norig, 1);
    typedef typename mp::make_homogeneous_tuple<T*, ntemps>::type
        temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, T* res,
            const T*& toright1, const T*& toright2)
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
template<class T, class Left, class right1_t, class right2_t, bool t1_ge_t2>
struct ternary_hhelper_1imm;
// Case where t1 >= t2
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper_1imm<T, Left, right1_t, right2_t, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef typename right1_t::ev_traits_t::temp_rule_t evr;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evr::temporaries_t::len;
    // t1 >= t2

    template<class Temps>
    static void doit(const Left& left, const right1_t& right1,
            Temps temps, T* res,
            const T*& toright1)
    {
        evl::doit(left._data(), mp::htuples::extract<t1>(temps), res);
        typename Temps::tail_t nores = mp::htuples::removeres(temps, res);
        evr::doit(right1._data(),
                mp::htuples::extract<t2>(nores), nores.head);
        toright1 = nores.head;
    }
};
// Case where t1 < t2
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper_1imm<T, Left, right1_t, right2_t, false>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef typename right1_t::ev_traits_t::temp_rule_t evr;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evr::temporaries_t::len;
    // t1 < t2

    template<class Temps>
    static void doit(const Left& left, const right1_t& right1,
            Temps temps, T* res,
            const T*& toright1)
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
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper<T, Left, right1_t, right2_t, false, true>
{
    typedef typename Left::ev_traits_t::temp_rule_t evl;
    typedef tools::evaluation_helper<right1_t> evhr1;
    static const unsigned t1 = evl::temporaries_t::len;
    static const unsigned t2 = evhr1::temporaries_t::len;
    static const unsigned ntemps = FLINT_MAX(2, FLINT_MAX(t1, t2) + (t1 == t2));

    typedef ternary_hhelper_1imm<T, Left, right1_t, right1_t, t1 >= t2> thh1;
    typedef typename mp::make_homogeneous_tuple<T*, ntemps>::type
        temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, T* res,
            const T*& toright1, const T*& toright2)
    {
        toright2 = &right2;
        thh1::doit(left, right1, temps, res, toright1);
    }
};

// Case where b is immediate.
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper<T, Left, right1_t, right2_t, true, false>
{
    typedef ternary_hhelper<T, Left, right2_t, right1_t, false, true> thh;
    typedef typename thh::temporaries_t temporaries_t;

    static void doit(const Left& left, const right1_t& right1,
            const right2_t& right2, temporaries_t temps, T* res,
            const T*& toright1, const T*& toright2)
    {
        thh::doit(left, right2, right1, temps, res, toright2, toright1);
    }
};

// Case where neither is immediate.
template<class T, class Left, class right1_t, class right2_t>
struct ternary_hhelper<T, Left, right1_t, right2_t, false, false>
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
    typedef typename mp::make_homogeneous_tuple<T*, ntemps>::type
        temporaries_t;

    struct resaccess
    {
        T* res;
        resaccess(T* r) : res(r) {};

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
        const T*& right;
        toaccess(const T*& r) : right(r) {};

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
            const right2_t& right2, temporaries_t temps, T* res,
            const T*& toright1, const T*& toright2)
    {
        // We re-order the temporaries in such a way that res is at the
        // very end. When evaluating things in the correct order, it is then
        // always correct to take temporaries from the front, and drop them
        // from the front.
        temporaries_t temps_reordered = mp::concat_tuple<
            typename temporaries_t::tail_t,
            typename mp::make_tuple<T*>::type>::doit(
                    mp::htuples::removeres(temps, res),
                    mp::make_tuple<T*>::make(res));
        doit_sort(left, right1, right2, temps_reordered,
                resaccess(res), toaccess(toright1), toaccess(toright2));
    }
};
} // tdetail
} // tools
} // flint

#endif
