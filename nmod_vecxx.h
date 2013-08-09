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

// TODO actual nmod_vecxx class
// TODO reference types
// TODO addmul

#ifndef NMOD_VECXX_H
#define NMOD_VECXX_H

#include <sstream>

#include "nmod_vec.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/stdmath.h"

namespace flint {
class nmodxx_ctx
{
private:
    nmod_t nmod;

public:
    const nmod_t& _nmod() const {return nmod;}
    explicit nmodxx_ctx(mp_limb_t n) {nmod_init(&nmod, n);}
    // no destruction necessary

    bool operator==(const nmodxx_ctx& o) const {return nmod.n == o.nmod.n;}
};

class nmodxx_ctx_srcref
{
private:
    const nmod_t& nmod;

    nmodxx_ctx_srcref(const nmod_t& nm) : nmod(nm) {}

public:
    const nmod_t& _nmod() const {return nmod;}
    nmodxx_ctx_srcref(const nmodxx_ctx& c) : nmod(c._nmod()) {}

    static nmodxx_ctx_srcref make(const nmod_t& nm)
        {return nmodxx_ctx_srcref(nm);}

    bool operator==(const nmodxx_ctx_srcref& o) const {return nmod.n == o.nmod.n;}
};

template<class Operation, class Data>
class nmodxx_expression
    : public expression<derived_wrapper<nmodxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::nmodxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS_NOFLINTCLASS(nmodxx_expression)
    FLINTXX_DEFINE_CTORS(nmodxx_expression)

    // static functions for nmodxx
    static nmodxx_expression make_nored(mp_limb_t n, nmodxx_ctx_srcref c)
    {
        return nmodxx_expression(Data::make_nored(n, c));
    }
    static nmodxx_expression red(mp_limb_t n, nmodxx_ctx_srcref c)
    {
        nmodxx_expression res = make_nored(n, c);
        res.reduce();
        return res;
    }
    // TODO more

    // only makes sense on immediates
    nmodxx_ctx_srcref _ctx() const {return this->_data().ctx;}
    const nmod_t& _nmod() const {return this->_data().ctx._nmod();}
    mp_limb_t& _limb() {return this->_data().limb;}
    const mp_limb_t& _limb() const {return this->_data().limb;}
    void reduce() {NMOD_RED(_limb(), _limb(), _nmod());}
    void set_nored(mp_limb_t n) {this->_data().limb = n;}

    nmodxx_ctx_srcref estimate_ctx() const;

    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx());
    }

    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
};

namespace detail {
struct nmodxx_data
{
    nmodxx_ctx_srcref ctx;
    mp_limb_t limb;

    nmodxx_data(nmodxx_ctx_srcref c) : ctx(c), limb(0) {}

private:
    nmodxx_data(mp_limb_t n, nmodxx_ctx_srcref c)
        : ctx(c), limb(n) {}

public:
    static nmodxx_data make_nored(mp_limb_t n, nmodxx_ctx_srcref c)
    {
        return nmodxx_data(n, c);
    }
};
} // detail

typedef nmodxx_expression<operations::immediate, detail::nmodxx_data> nmodxx;

#if 0
namespace traits {
// Temporary merging isn't really any use here
template<> struct use_temporary_merging<nmodxx> : mp::false_ { };
} // traits
#endif

namespace detail {
struct is_nmodxx_predicate
{
    // TODO this needs change with reference types
    template<class T> struct type : mp::equal_types<T, nmodxx> { };
};
} // detail
template<class Operation, class Data>
inline nmodxx_ctx_srcref
nmodxx_expression<Operation, Data>::estimate_ctx() const
{
    return tools::find_subexpr<detail::is_nmodxx_predicate>(*this)._ctx();
}

namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct NMODXX_COND_S : mp::equal_types<T, nmodxx> { };
#define NMODXX_COND_T NMODXX_COND_S

// TODO references
FLINT_DEFINE_GET(equals, bool, nmodxx, e1._limb() == e2._limb())

FLINT_DEFINE_GET_COND(conversion, mp_limb_t, NMODXX_COND_S, from._limb())

template<class Nmod>
struct to_string<Nmod, typename mp::enable_if< NMODXX_COND_S<Nmod> >::type>
{
    static std::string get(const nmodxx& i, int base /* ignored */)
    {
        std::ostringstream oss;
        oss << i._limb() << " mod " << i._nmod().n;
        return oss.str();
    }
};

FLINT_DEFINE_DOIT_COND2(assignment, NMODXX_COND_T, NMODXX_COND_S,
        to._limb() = from._limb())

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, nmodxx, NMODXX_COND_S, NMODXX_COND_S,
        to.set_nored(nmod_add(e1._limb(), e2._limb(), to._nmod())))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, nmodxx, NMODXX_COND_S, NMODXX_COND_S,
        to.set_nored(nmod_mul(e1._limb(), e2._limb(), to._nmod())))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, nmodxx, NMODXX_COND_S, NMODXX_COND_S,
        to.set_nored(nmod_sub(e1._limb(), e2._limb(), to._nmod())))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, nmodxx, NMODXX_COND_S, NMODXX_COND_S,
        to.set_nored(nmod_div(e1._limb(), e2._limb(), to._nmod())))

FLINT_DEFINE_UNARY_EXPR_COND(negate, nmodxx, NMODXX_COND_S,
        to.set_nored(nmod_neg(from._limb(), to._nmod())))

FLINT_DEFINE_UNARY_EXPR_COND(inv_op, nmodxx, NMODXX_COND_S,
        to.set_nored(nmod_inv(from._limb(), to._nmod())))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, nmodxx, NMODXX_COND_S,
        traits::is_unsigned_integer,
        to.set_nored(nmod_pow_ui(e1._limb(), e2, to._nmod())))
}
} // flint

#endif
