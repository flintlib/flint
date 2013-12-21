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

#ifndef QADICXX_H
#define QADICXX_H

#include <algorithm> // std::max

#include "qadic.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/matrix.h" // trace ...

#include "padicxx.h"

namespace flint {
FLINT_DEFINE_BINOP(frobenius)
FLINT_DEFINE_UNOP(norm)
FLINT_DEFINE_UNOP(norm_analytic)
FLINT_DEFINE_UNOP(norm_resultant)

class qadicxx_ctx
{
private:
    mutable qadic_ctx_t ctx;

public:
    // NB: you must not modify user-visible state of ctx through a constant
    // instance of qadicxx_ctx
    qadic_ctx_t& _ctx() const {return ctx;}
    padicxx_ctx_srcref pctx() const
        {return padicxx_ctx_srcref::make(&ctx->pctx);}

    // TODO more constructors? Should we wrap padic_print_mode?
    qadicxx_ctx(fmpzxx_srcref p, slong d, slong min, slong max,
        padic_print_mode mode, const char* var = "x")
    {
        qadic_ctx_init_conway(ctx, p._fmpz(), d, min, max, var, mode);
    }

    ~qadicxx_ctx() {qadic_ctx_clear(ctx);}
};

inline void print(const qadicxx_ctx& c)
{
    qadic_ctx_print(c._ctx());
}

namespace traits {
template<class T> struct has_qadicxx_ctx : mp::false_ { };

template<class T> struct is_qadic_expr
    : has_qadicxx_ctx<typename tools::evaluation_helper<T>::type> { };
} // traits

namespace detail {
struct has_qadicxx_ctx_predicate
{
    template<class T> struct type : traits::has_qadicxx_ctx<T> { };
};
} // detail

namespace tools {
template<class Expr>
const qadicxx_ctx& find_qadicxx_ctx(const Expr& e)
{
    return tools::find_subexpr<detail::has_qadicxx_ctx_predicate>(e).get_qctx();
}
} // tools

namespace detail {
template<class Qadic>
struct qadic_traits
{
    static slong prec(const Qadic& q) {return tools::padic_output_prec(q);}
};
} // detail

template<class Operation, class Data>
class qadicxx_expression
    : public expression<derived_wrapper<qadicxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::qadicxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(qadicxx_expression)
    FLINTXX_DEFINE_CTORS(qadicxx_expression)
    FLINTXX_DEFINE_C_REF(qadicxx_expression, qadic_struct, _qadic)

    typedef detail::qadic_traits<qadicxx_expression> traits_t;

    // These only make sense with immediates
    void reduce() {qadic_reduce(_qadic(), _ctx());}
    void set_zero() {qadic_zero(_qadic());}
    void set_one() {qadic_one(_qadic());}
    void set_gen() {qadic_gen(_qadic(), _ctx());}

    const qadicxx_ctx& get_qctx() const {return this->_data().ctx;}
    padicxx_ctx_srcref get_ctx() const {return get_qctx().pctx();}
    qadic_ctx_t& _ctx() const {return get_qctx()._ctx();}

    // these only make sense with qadicxx
    static qadicxx_expression zero(const qadicxx_ctx& ctx)
        {return qadicxx_expression(ctx);}
    static qadicxx_expression zero(const qadicxx_ctx& ctx, slong N)
        {return qadicxx_expression(ctx, N);}
    static qadicxx_expression one(const qadicxx_ctx& ctx)
    {
        qadicxx_expression res(ctx);
        res.set_one();
        return res;
    }
    static qadicxx_expression one(const qadicxx_ctx& ctx, slong N)
    {
        qadicxx_expression res(ctx, N);
        res.set_one();
        return res;
    }
    static qadicxx_expression gen(const qadicxx_ctx& ctx)
    {
        qadicxx_expression res(ctx);
        res.set_gen();
        return res;
    }
    static qadicxx_expression gen(const qadicxx_ctx& ctx, slong N)
    {
        qadicxx_expression res(ctx, N);
        res.set_gen();
        return res;
    }

    template<class Padic>
    static qadicxx_expression from_ground(const qadicxx_ctx& ctx,
            const Padic& p,
            typename mp::enable_if<traits::is_padicxx<Padic> >::type* = 0)
    {
        qadicxx_expression res(ctx);
        res = p;
        return res;
    }

    template<class Padic>
    static qadicxx_expression from_ground(const qadicxx_ctx& ctx, slong N,
            const Padic& p,
            typename mp::enable_if<traits::is_padicxx<Padic> >::type* = 0)
    {
        qadicxx_expression res(ctx, N);
        res = p;
        return res;
    }

    static qadicxx_expression randtest(frandxx& state,
            const qadicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        qadicxx_expression res(ctx, prec);
        qadic_randtest(res._qadic(), state._data(), ctx._ctx());
        return res;
    }
    static qadicxx_expression randtest_not_zero(frandxx& state,
            const qadicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        qadicxx_expression res(ctx, prec);
        qadic_randtest_not_zero(res._qadic(), state._data(), ctx._ctx());
        return res;
    }
    static qadicxx_expression randtest_val(frandxx& state,
            slong val, const qadicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        qadicxx_expression res(ctx, prec);
        qadic_randtest_val(res._qadic(), state._data(), val, ctx._ctx());
        return res;
    }
    static qadicxx_expression randtest_int(frandxx& state,
            const qadicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        qadicxx_expression res(ctx, prec);
        qadic_randtest_int(res._qadic(), state._data(), ctx._ctx());
        return res;
    }

    const qadicxx_ctx& estimate_ctx() const
    {
        return tools::find_qadicxx_ctx(*this);
    }

    // Create a temporary. The context will be estimated, and the precision
    // will be the maximum of all subexpressions.
    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx(), prec());
    }

    // TODO randomisation

    typename flint_classes::to_srcref<typename base_t::derived_t>::type
    toN(slong N) const
    {
        return flint_classes::to_srcref<typename base_t::derived_t>::type::make(
                this->_data().inner, get_qctx(), N);
    }

    slong prec () const {return traits_t::prec(*this);}

    // these cause evaluation
    slong val() const {return qadic_val(this->evaluate()._qadic());}
    bool is_zero() const {return qadic_is_zero(this->evaluate()._qadic());}
    bool is_one() const {return qadic_is_one(this->evaluate()._qadic());}

    // forwarding of lazy functions
    FLINTXX_DEFINE_MEMBER_BINOP(frobenius)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_UNOP(exp)
    FLINTXX_DEFINE_MEMBER_UNOP(exp_balanced)
    FLINTXX_DEFINE_MEMBER_UNOP(exp_rectangular)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(log)
    FLINTXX_DEFINE_MEMBER_UNOP(log_balanced)
    FLINTXX_DEFINE_MEMBER_UNOP(teichmuller)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(padicxx, trace)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(padicxx, norm)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(padicxx, norm_analytic)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(padicxx, norm_resultant)
};

namespace detail {
struct qadic_data;
}

typedef qadicxx_expression<operations::immediate, detail::qadic_data> qadicxx;
typedef qadicxx_expression<operations::immediate,
            flint_classes::ref_data<qadicxx, qadic_struct> > qadicxx_ref;
typedef qadicxx_expression<operations::immediate, flint_classes::srcref_data<
    qadicxx, qadicxx_ref, qadic_struct> > qadicxx_srcref;

namespace traits {
template<> struct has_padicxx_ctx<qadicxx> : mp::true_ { };
template<> struct has_padicxx_ctx<qadicxx_ref> : mp::true_ { };
template<> struct has_padicxx_ctx<qadicxx_srcref> : mp::true_ { };
template<> struct has_qadicxx_ctx<qadicxx> : mp::true_ { };
template<> struct has_qadicxx_ctx<qadicxx_ref> : mp::true_ { };
template<> struct has_qadicxx_ctx<qadicxx_srcref> : mp::true_ { };

template<class T> struct is_qadicxx : flint_classes::is_Base<padicxx, T> { };
} // traits

namespace detail {
template<>
struct qadic_traits<qadicxx_srcref>
{
    template<class Q>
    static slong prec(const Q& q) {return q._data().N;}
};
template<>
struct qadic_traits<qadicxx_ref>
{
    template<class Q>
    static slong prec(const Q& q) {return qadic_prec(q._qadic());}
};
template<> struct qadic_traits<qadicxx> : qadic_traits<qadicxx_ref> { };
}

PADICXX_DEFINE_REF_STRUCTS_(qadicxx, qadic_struct, qadic_prec, const qadicxx_ctx&)

namespace detail {
struct qadic_data
{
    typedef qadic_t& data_ref_t;
    typedef const qadic_t& data_srcref_t;

    const qadicxx_ctx& ctx;
    qadic_t inner;

    qadic_data(const qadicxx_ctx& c)
        : ctx(c)
    {
        qadic_init(inner);
    }

    qadic_data(const qadicxx_ctx& c, slong N)
        : ctx(c)
    {
        qadic_init2(inner, N);
    }

    qadic_data(const qadic_data& o)
        : ctx(o.ctx)
    {
        qadic_init2(inner, qadic_prec(o.inner));
        qadic_set(inner, o.inner, ctx._ctx());
    }

    ~qadic_data() {qadic_clear(inner);}

    qadic_data(qadicxx_srcref c)
        : ctx(c.get_qctx())
    {
        qadic_init2(inner, c.prec());
        qadic_set(inner, c._qadic(), ctx._ctx());
    }
};
} // detail

#define QADICXX_COND_S FLINTXX_COND_S(qadicxx)
#define QADICXX_COND_T FLINTXX_COND_T(qadicxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, QADICXX_COND_T, QADICXX_COND_S,
        qadic_set(to._qadic(), from._qadic(), to._ctx()))

FLINT_DEFINE_DOIT_COND2(assignment, QADICXX_COND_T, traits::is_unsigned_integer, 
        qadic_set_ui(to._qadic(), from, to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, QADICXX_COND_T, PADICXX_COND_S,
        padic_poly_set_padic(to._qadic(), from._padic(), from._ctx()))

FLINT_DEFINE_PRINT_PRETTY_COND(QADICXX_COND_S,
        qadic_fprint_pretty(to, from._qadic(), from._ctx()))

template<class T>
struct conversion<padicxx, T,
    typename mp::enable_if< QADICXX_COND_S<T> >::type>
{
    static padicxx get(const T& from)
    {
        padicxx res(from.estimate_ctx().pctx(), from.prec());
        execution_check(qadic_get_padic(res._padic(), from._qadic(), from._ctx()),
                "get_padic", "qadic");
        return res;
    }
};

FLINTXX_DEFINE_SWAP(qadicxx, qadic_swap(e1._qadic(), e2._qadic()))

FLINTXX_DEFINE_EQUALS(qadicxx, qadic_equal(e1._qadic(), e2._qadic()))

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, qadicxx, QADICXX_COND_S, QADICXX_COND_S,
        qadic_add(to._qadic(), e1._qadic(), e2._qadic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, qadicxx, QADICXX_COND_S, QADICXX_COND_S,
        qadic_sub(to._qadic(), e1._qadic(), e2._qadic(), to._ctx()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, qadicxx, QADICXX_COND_S, QADICXX_COND_S,
        qadic_mul(to._qadic(), e1._qadic(), e2._qadic(), to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, qadicxx, QADICXX_COND_S,
        qadic_neg(to._qadic(), from._qadic(), to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(inv_op, qadicxx, QADICXX_COND_S,
        qadic_inv(to._qadic(), from._qadic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, qadicxx, QADICXX_COND_S,
        traits::is_fmpzxx,
        qadic_pow(to._qadic(), e1._qadic(), e2._fmpz(), to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(exp_op, qadicxx, QADICXX_COND_S,
        execution_check(
            qadic_exp(to._qadic(), from._qadic(), to._ctx()), "exp", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(exp_balanced_op, qadicxx, QADICXX_COND_S,
        execution_check(qadic_exp_balanced(
                to._qadic(), from._qadic(), to._ctx()), "exp_balanced", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(exp_rectangular_op, qadicxx, QADICXX_COND_S,
        execution_check(qadic_exp_rectangular(
                to._qadic(), from._qadic(), to._ctx()),
            "exp_rectangular", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_op, qadicxx, QADICXX_COND_S,
        execution_check(
            qadic_log(to._qadic(), from._qadic(), to._ctx()), "log", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_rectangular_op, qadicxx, QADICXX_COND_S,
        execution_check(qadic_log_rectangular(
                to._qadic(), from._qadic(), to._ctx()),
            "log_rectangular", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_balanced_op, qadicxx, QADICXX_COND_S,
        execution_check(qadic_log_balanced(
                to._qadic(), from._qadic(), to._ctx()), "log_balanced", "qadic"))
FLINT_DEFINE_UNARY_EXPR_COND(teichmuller_op, qadicxx, QADICXX_COND_S,
            qadic_teichmuller(to._qadic(), from._qadic(), to._ctx()))

FLINT_DEFINE_BINARY_EXPR_COND2(frobenius_op, qadicxx, QADICXX_COND_S,
        traits::fits_into_slong,
        qadic_frobenius(to._qadic(), e1._qadic(), e2, to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(trace_op, padicxx, QADICXX_COND_S,
            qadic_trace(to._padic(), from._qadic(), from._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(norm_op, padicxx, QADICXX_COND_S,
            qadic_norm(to._padic(), from._qadic(), from._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(norm_analytic_op, padicxx, QADICXX_COND_S,
            qadic_norm_analytic(to._padic(), from._qadic(), from._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(norm_resultant_op, padicxx, QADICXX_COND_S,
            qadic_norm_resultant(to._padic(), from._qadic(), from._ctx()))
} // rules
} // flint

#endif
