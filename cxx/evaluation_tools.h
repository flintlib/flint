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
template<class T, class U>
struct symmetric_cmp
    : mp::if_<traits::is_implemented<rules::cmp<T, U> >,
          rules::cmp<T, U>,
          typename mp::if_<traits::is_implemented<rules::cmp<U, T> >,
              tdetail::cmp_invert<T, U>,
              rules::UNIMPLEMENTED
            >::type
        >::type { };

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
    typedef typename Expr1::ev_traits_t::rule_t rule1_t;
    typedef typename Expr2::ev_traits_t::rule_t rule2_t;

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

    // TODO sometimes the other way round uses fewer temporaries
    typedef mp::merge_tuple<typename mp::make_tuple<return2_t*>::type,
                 temporaries2_t> merger2;
    typedef mp::merge_tuple<tuple<return1_t*, typename merger2::type>,
                 temporaries1_t> merger1;

    return1_t* ret1;
    return2_t* ret2;

public:
    typedef typename merger1::type temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, const Expr2& e2)
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

    ref1_t get1() {return *ret1;}
    ref2_t get2() {return *ret2;}
};

///////////////////////////////////////////////////////////////////////////
// Helper to evaluate three homogeneous terms
///////////////////////////////////////////////////////////////////////////
} // tools
} // flint

#endif
