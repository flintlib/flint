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

#include <iostream>

#include "../flint.h" // FLINT_MAX and FLINT_MIN

#include "expression_traits.h"
#include "mp.h"
#include "rules.h"
#include "tuple.h"

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

// A version of equals which uses cmp if possible
namespace tdetail {
template<class T, class U, class Enable = void>
struct equals_using_cmp_ : rules::UNIMPLEMENTED { };
template<class T, class U>
struct equals_using_cmp_<T, U,
    typename mp::enable_if<
        traits::is_implemented<symmetric_cmp<T, U> > >::type>
{
    static bool get(const T& t, const U& u)
    {
        return tools::symmetric_cmp<T, U>::get(t, u) == 0;
    }
};
} // tdetail
template<class T, class U>
struct equals_using_cmp
    : mp::if_<traits::is_implemented<rules::equals<T, U> >,
        rules::equals<T, U>, tdetail::equals_using_cmp_<T, U> >::type { };

// Automatic printing if to_string is implemented
namespace tdetail {
template<class T, class Enable = void>
struct print_using_str_ : rules::UNIMPLEMENTED { };
template<class T>
struct print_using_str_<T,
    typename mp::enable_if<traits::is_implemented<rules::to_string<T> > >::type>
{
    static void doit(const T& v, std::ostream& o)
    {
        int base = 10;
        std::ios_base::fmtflags ff = o.flags();
        if(ff & o.hex)
            base = 16;
        if(ff & o.oct)
            base = 8;
        o << v.to_string(base);
    }
};
} // tdetail
template<class T>
struct print_using_str
    : mp::if_<traits::is_implemented<rules::print<T> >,
        rules::print<T>, tdetail::print_using_str_<T> >::type { };

// Finding a subexpression of precsribed type
namespace tdetail {
template<class Pred, class Data, class Enable = void>
struct find_subexpr_helper2;

template<class Pred, class Expr, class Enable = void>
struct find_subexpr_helper
{
    typedef find_subexpr_helper2<Pred, typename Expr::data_t> fsh;
    typedef typename fsh::rtype rtype;
    static const bool val = fsh::val;
    static rtype get(const Expr& e) {return fsh::get(e._data());}
};

template<class Pred, class Expr>
struct find_subexpr_helper<Pred, Expr,
    typename mp::enable_if<typename Pred::template type<Expr> >::type>
{
    static const bool val = true;
    typedef const Expr& rtype;
    static rtype get(rtype t) {return t;}
};

template<class Pred, class Expr>
struct find_subexpr_helper<Pred, Expr,
    typename mp::enable_if<mp::and_<traits::is_immediate<Expr>,
        mp::not_<typename Pred::template type<Expr> > > >::type>
{
    static const bool val = false;
    typedef void rtype;
};

template<class Pred, class Data, class Enable>
struct find_subexpr_helper2
{
    typedef find_subexpr_helper2<Pred, typename Data::tail_t> fsh;
    typedef typename fsh::rtype rtype;
    static const bool val = fsh::val;
    static rtype get(const Data& d) {return fsh::get(d.tail);}
};

template<class Pred, class Head, class Tail>
struct find_subexpr_helper2<Pred, tuple<Head, Tail>,
    typename mp::enable_if<find_subexpr_helper<Pred,
        typename traits::basetype<Head>::type> >::type>
{
    static const bool val = true;
    typedef typename traits::basetype<Head>::type head_t;
    typedef find_subexpr_helper<Pred, head_t> fsh;
    typedef typename fsh::rtype rtype;
    static rtype get(const tuple<Head, Tail>& d) {return fsh::get(d.head);}
};

template<class Pred>
struct find_subexpr_helper2<Pred, empty_tuple>
{
    static const bool val = false;
    typedef void rtype;
};
} // tdetail

// A predicate which applies if the argument type equals T.
template<class T>
struct equal_types_pred
{
    template<class U> struct type : mp::equal_types<T, U> { };
};

// Given an expression template Expr, traverse the tree of data arguments
// until an argument matching the predicate Pred is found. Here pred must have
// a member template "type" performing the boolean computation.
// See equal_types_pred for an example.
// If there is no matching subexpression, a compile time error will be
// encountered.
// The current implementation performs depth-first search.
template<class Pred, class Expr>
inline typename tdetail::find_subexpr_helper<Pred, Expr>::rtype
find_subexpr(const Expr& e)
{
    return tdetail::find_subexpr_helper<Pred, Expr>::get(e);
}

// Find a subexpression of type T.
template<class T, class Expr>
inline const T& find_subexpr_T(const Expr& e)
{
    return find_subexpr<equal_types_pred<T> >(e);
}

// Boolean computation to determine if find_subexpr above will work.
template<class Pred, class Expr>
struct has_subexpr
    : tdetail::find_subexpr_helper<Pred, Expr> { };


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
    typedef ftype etype;
    static ftype get(const type& t) {return t;}

    typedef empty_tuple temporaries_t;
};

template<class T>
struct evaluation_helper<T,
    typename mp::enable_if<traits::is_lazy_expr<T> >::type>
{
    typedef typename T::evaluated_t type;
    typedef const typename T::evaluated_t& ftype;
    typedef type etype;
    static type get(const T& t) {return t.evaluate();}

    typedef typename T::ev_traits_t::temp_rule_t::temporaries_t temporaries_t;
};


///////////////////////////////////////////////////////////////////////////
// Helper to evaluate n terms
///////////////////////////////////////////////////////////////////////////
// The template argument is an arbitrary argument tuple
template<class Args, class Enable = void>
struct evaluate_n;

// Count the number of non-immediate terms in arguments
template<class Args>
struct count_nonimm
{
    static const unsigned val =
        traits::is_lazy_expr<typename Args::head_t>::val
        + count_nonimm<typename Args::tail_t>::val;
};
template<>
struct count_nonimm<empty_tuple>
{
    static const unsigned val = 0;
};

template<class Args>
struct evaluated_args_tuple
{
    typedef typename tools::evaluation_helper<typename Args::head_t>::ftype evt;
    typedef typename evaluated_args_tuple<typename Args::tail_t>::type tail_t;
    typedef tuple<evt, tail_t> type;
};
template<>
struct evaluated_args_tuple<empty_tuple>
{
    typedef empty_tuple type;
};

namespace tdetail {
// Unoptimized evaluation (in order)
// The constructor sets up any local temporaries (coming from non-merging).
// Then init does the actual computation. Note that init chains at the end
// (in contrast to construction), so we do evaluate in order.
template<class Args, class Enable = void>
struct evaluate_n_unopt
{
    // case where head is immediate
    typedef evaluate_n_unopt<typename Args::tail_t> next_t;

    typedef typename next_t::temporaries_t temporaries_t;
    typedef typename Args::head_t headr_t;

    next_t next;
    headr_t res;

    headr_t gethead() const
    {
        return res;
    }

    void init(const Args& args, temporaries_t temps)
    {
        next.init(args.tail, temps);
    }
    evaluate_n_unopt(const Args& args)
        : next(args.tail), res(args.head) {}
};
template<>
struct evaluate_n_unopt<empty_tuple>
{
    // basecase
    typedef empty_tuple temporaries_t;

    void init(empty_tuple, empty_tuple) {}
    evaluate_n_unopt(empty_tuple) {}
};
template<class Args>
struct evaluate_n_unopt<Args, typename mp::enable_if<mp::and_<
    traits::is_lazy_expr<typename Args::head_t>,
    mp::not_<traits::use_temporary_merging<
        typename Args::head_t::evaluated_t> > > >::type>
{
    // Case with non-merging lazy head
    typedef evaluate_n_unopt<typename Args::tail_t> next_t;
    typedef typename Args::head_t expr_t;
    typedef typename expr_t::ev_traits_t::temp_rule_t rule_t;
    typedef typename rule_t::return_t tmp_t;
    typedef mp::merge_tuple<typename rule_t::temporaries_t,
              typename next_t::temporaries_t> merger;

    typedef typename merger::type temporaries_t;
    typedef typename traits::forwarding<tmp_t>::type headr_t;

    next_t next;
    tmp_t tmp;

    headr_t gethead() const
    {
        return tmp;
    }

    evaluate_n_unopt(const Args& args)
        : next(args.tail), tmp(
                rules::instantiate_temporaries<expr_t, tmp_t>::get(args.head))
    {}
    void init(const Args& args, temporaries_t temps)
    {
        rule_t::doit(args.head._data(), merger::get_first(temps), &tmp);
        next.init(args.tail, merger::get_second(temps));
    }
};
template<class Args>
struct evaluate_n_unopt<Args, typename mp::enable_if<mp::and_<
    traits::is_lazy_expr<typename Args::head_t>,
    traits::use_temporary_merging<
        typename Args::head_t::evaluated_t> > >::type>
{
    // Case with merging lazy head
    typedef evaluate_n_unopt<typename Args::tail_t> next_t;
    typedef typename Args::head_t expr_t;
    typedef typename expr_t::ev_traits_t::temp_rule_t rule_t;
    typedef typename rule_t::return_t tmp_t;
    typedef mp::merge_tuple<typename rule_t::temporaries_t,
              tuple<tmp_t*, typename next_t::temporaries_t> > merger;

    typedef typename merger::type temporaries_t;
    typedef typename traits::forwarding<tmp_t>::type headr_t;

    next_t next;
    tmp_t* tmp;

    headr_t gethead() const
    {
        return *tmp;
    }

    evaluate_n_unopt(const Args& args)
        : next(args.tail) {}
    void init(const Args& args, temporaries_t temps)
    {
        tmp = merger::get_second(temps).head;
        rule_t::doit(args.head._data(), merger::get_first(temps), tmp);
        next.init(args.tail, merger::get_second(temps).tail);
    }
};

template<class Args, unsigned n>
struct unopt_get
{
    typedef unopt_get<typename Args::tail_t, n-1> getn;
    typedef typename getn::type type;
    static type get(const evaluate_n_unopt<Args>& e) {return getn::get(e.next);}
};
template<class Args>
struct unopt_get<Args, 0>
{
    typedef evaluate_n_unopt<Args> evalt;
    typedef typename evalt::headr_t type;
    static type get(const evalt& e) {return e.gethead();}
};

template<class Args>
struct unopt_gettuple
{
    typedef unopt_gettuple<typename Args::tail_t> next;
    typedef evaluate_n_unopt<Args> eval_t;

    typedef tuple<typename eval_t::headr_t, typename next::type> type;
    static type get(const eval_t& e)
    {
        return type(e.gethead(), next::get(e.next));
    }
};
template<>
struct unopt_gettuple<empty_tuple>
{
    typedef empty_tuple type;
    template<class T>
    static type get(const T&) {return empty_tuple();}
};

// Optimized case with precisely two non-immediates
template<class Args, unsigned found = 0, class Enable = void>
struct evaluate_n_2_analyze
{
    typedef evaluate_n_2_analyze<typename Args::tail_t, found> next;
    static const unsigned first = next::first + 1;
    static const unsigned second = next::second + 1;
};
template<class Args>
struct evaluate_n_2_analyze<Args, 0,
    typename mp::enable_if<traits::is_lazy_expr<typename Args::head_t> >::type>
{
    typedef evaluate_n_2_analyze<typename Args::tail_t, 1> next;
    static const unsigned first = 0;
    static const unsigned second = next::second + 1;
};
template<class Args>
struct evaluate_n_2_analyze<Args, 1,
    typename mp::enable_if<traits::is_lazy_expr<typename Args::head_t> >::type>
{
    static const unsigned first = 0;
    static const unsigned second = 0;
};

template<class Expr1, class Expr2, class Enable = void>
struct evaluate_2;
// Case where neither is immediate, no merging
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_lazy_expr<Expr2>,
        mp::not_<traits::use_temporary_merging<typename Expr1::evaluated_t> >,
        mp::not_<traits::use_temporary_merging<typename Expr2::evaluated_t> >
      > >::type>
{
private:
    typedef typename Expr1::ev_traits_t::temp_rule_t rule1_t;
    typedef typename Expr2::ev_traits_t::temp_rule_t rule2_t;

public:
    typedef typename rule1_t::return_t return1_t;
    typedef typename rule2_t::return_t return2_t;

private:
    typedef typename rule1_t::temporaries_t temporaries1_t;
    typedef typename rule2_t::temporaries_t temporaries2_t;

    typedef mp::merge_tuple<temporaries1_t, temporaries2_t> merger;
    return1_t tmp1;
    return2_t tmp2;

public:
    typedef typename merger::type temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, const Expr2& e2)
        : tmp1(rules::instantiate_temporaries<Expr1, return1_t>::get(e1)),
          tmp2(rules::instantiate_temporaries<Expr2, return2_t>::get(e2))
    {
        rule1_t::doit(e1._data(), merger::get_first(temps), &tmp1);
        rule2_t::doit(e2._data(), merger::get_second(temps), &tmp2);
    }

    const return1_t& get1() const {return tmp1;}
    const return2_t& get2() const {return tmp2;}
};

// Case where neither is immediate, first has merging, second does not
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_lazy_expr<Expr2>,
        traits::use_temporary_merging<typename Expr1::evaluated_t>,
        mp::not_<traits::use_temporary_merging<typename Expr2::evaluated_t> >
      > >::type>
{
private:
    typedef typename Expr1::ev_traits_t::temp_rule_t rule1_t;
    typedef typename Expr2::ev_traits_t::temp_rule_t rule2_t;

public:
    typedef typename rule1_t::return_t return1_t;
    typedef typename rule2_t::return_t return2_t;

private:
    typedef typename rule1_t::temporaries_t temporaries1_t;
    typedef typename rule2_t::temporaries_t temporaries2_t;

    typedef mp::merge_tuple<typename mp::make_tuple<return1_t*>::type,
              temporaries1_t> merger1;
    typedef mp::merge_tuple<typename merger1::type, temporaries2_t> merger2;
    return2_t tmp2;
    return1_t* ret1;
    
public:
    typedef typename merger2::type temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, const Expr2& e2)
         : tmp2(rules::instantiate_temporaries<Expr2, return2_t>::get(e2))
    {
        rule2_t::doit(e2._data(), merger2::get_second(temps), &tmp2);
        ret1 = merger1::get_first(merger2::get_first(temps)).head;
        rule1_t::doit(e1._data(),
                merger1::get_second(merger2::get_first(temps)), ret1);
    }

    const return1_t& get1() const {return *ret1;}
    const return2_t& get2() const {return tmp2;}
};

// Case where neither is immediate, second has merging, first does not
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_lazy_expr<Expr2>,
        traits::use_temporary_merging<typename Expr2::evaluated_t>,
        mp::not_<traits::use_temporary_merging<typename Expr1::evaluated_t> >
      > >::type>
{
    // XXX this is copy-paste from above case where right is immediate
private:
    typedef evaluate_2<Expr2, Expr1> ev2_t;
    ev2_t ev2;

public:
    typedef typename ev2_t::return1_t return2_t;
    typedef typename ev2_t::return2_t return1_t;
    typedef typename ev2_t::temporaries_t temporaries_t;

    evaluate_2(temporaries_t temps, const Expr1& e1, const Expr2& e2)
        : ev2(temps, e2, e1) {};
    const return1_t& get1() const {return ev2.get2();}
    const return2_t& get2() const {return ev2.get1();}
};

// Case where neither is immediate, all merging
template<class Expr1, class Expr2>
struct evaluate_2<Expr1, Expr2,
    typename mp::enable_if<mp::and_<
        traits::is_lazy_expr<Expr1>,
        traits::is_lazy_expr<Expr2>,
        traits::use_temporary_merging<typename Expr1::evaluated_t>,
        traits::use_temporary_merging<typename Expr2::evaluated_t> > >::type>
{
private:
    typedef typename Expr1::ev_traits_t::temp_rule_t rule1_t;
    typedef typename Expr2::ev_traits_t::temp_rule_t rule2_t;

public:
    typedef typename rule1_t::return_t return1_t;
    typedef typename rule2_t::return_t return2_t;

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

    const return1_t& get1() const {return *ret1;}
    const return2_t& get2() const {return *ret2;}
};

template<unsigned first, unsigned second, unsigned n>
struct evaluate_n_2_get
{
    template<class Args, class First, class Second>
    static typename mp::tuple_get<Args, n>::type get(const Args& args,
            const First&, const Second&)
    {
        return mp::tuple_get<Args, n>::get(args);
    }
};
template<unsigned first, unsigned second>
struct evaluate_n_2_get<first, second, first>
{
    template<class Args, class First, class Second>
    static const First& get(const Args&, const First& f, const Second&)
    {
        return f;
    }
};
template<unsigned first, unsigned second>
struct evaluate_n_2_get<first, second, second>
{
    template<class Args, class First, class Second>
    static const Second& get(const Args&, const First&, const Second& s)
    {
        return s;
    }
};

template<class Tuple, unsigned n = 0>
struct evaluate_n_2_gettuple
{
    template<class T>
    static Tuple get(const T& t)
    {
        return Tuple(t.template get<n>(),
                evaluate_n_2_gettuple<typename Tuple::tail_t, n+1>::get(t));
    }
};
template<unsigned n>
struct evaluate_n_2_gettuple<empty_tuple, n>
{
    template<class T>
    static empty_tuple get(const T&) {return empty_tuple();}
};
} // tdetail

template<class Args, class Enable>
struct evaluate_n
{
    typedef tdetail::evaluate_n_unopt<Args> eval_t;
    typedef typename eval_t::temporaries_t temporaries_t;
    typedef typename tdetail::unopt_gettuple<Args>::type evtup_t;

    eval_t eval;

    evaluate_n(const Args& args, temporaries_t temps)
        : eval(args)
    {
        eval.init(args, temps);
    }

    template<unsigned n>
    typename tdetail::unopt_get<Args, n>::type get() const
    {
        return tdetail::unopt_get<Args, n>::get(eval);
    }

    evtup_t gettuple() const
    {
        return tdetail::unopt_gettuple<Args>::get(eval);
    }
};

template<class Args>
struct evaluate_n<Args,
    typename mp::enable_if_v<count_nonimm<Args>::val == 2>::type>
{
    typedef tdetail::evaluate_n_2_analyze<Args> analysis;
    static const unsigned first = analysis::first;
    static const unsigned second = analysis::second;

    typedef mp::tuple_get<Args, first> getfirst;
    typedef mp::tuple_get<Args, second> getsecond;
    typedef typename getfirst::type first_t;
    typedef typename getsecond::type second_t;

    typedef tdetail::evaluate_2<first_t, second_t> ev_t;
    typedef typename ev_t::temporaries_t temporaries_t;
    const Args& args;
    ev_t ev;

    evaluate_n(const Args& a, temporaries_t temps)
        : args(a), ev(temps, getfirst::get(a), getsecond::get(a)) {}

    typedef typename evaluated_args_tuple<Args>::type evtup_t;

    template<unsigned n>
    typename mp::tuple_get<evtup_t, n>::type get() const
    {
        return tdetail::evaluate_n_2_get<first, second, n>::get(
                args, ev.get1(), ev.get2());
    }

    evtup_t gettuple() const
    {
        return tdetail::evaluate_n_2_gettuple<evtup_t>::get(*this);
    }
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
// NOTE: This implementation does not honor use_temporary_merging!
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

// A helper condition for use with FLINT_DEFINE_*_COND?
template<class T>
struct is_bool : mp::equal_types<T, bool> { };
} // tools
} // flint

#endif
