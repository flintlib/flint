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

// TODO reference types
// TODO addmul

// TODO document nmod_vecxx

#ifndef NMOD_VECXX_H
#define NMOD_VECXX_H

#include <sstream>

#include "nmod_vec.h"

// TODO reduce dependencies?
#include "fmpzxx.h"
#include "fmpqxx.h"

#include "flintxx/expression.h"
#include "flintxx/evaluation_tools.h"
#include "flintxx/flint_classes.h"
#include "flintxx/stdmath.h"
#include "flintxx/vector.h"

namespace flint {
//////////////////////////////////////////////////////////////////////////////
// NMOD CLASS AND RULES
//////////////////////////////////////////////////////////////////////////////
class nmodxx_ctx
{
private:
    nmod_t nmod;

public:
    const nmod_t& _nmod() const {return nmod;}
    explicit nmodxx_ctx(mp_limb_t n) {nmod_init(&nmod, n);}
    // no destruction necessary

    bool operator==(const nmodxx_ctx& o) const {return nmod.n == o.nmod.n;}
    mp_limb_t n() const {return nmod.n;}
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
    mp_limb_t n() const {return nmod.n;}
};

namespace detail {
struct nmodxx_fake_c_type { };
} // detail
template<class Operation, class Data>
class nmodxx_expression
    : public expression<derived_wrapper<nmodxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::nmodxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(nmodxx_expression)
    FLINTXX_DEFINE_CTORS(nmodxx_expression)
    FLINTXX_DEFINE_C_REF(nmodxx_expression, detail::nmodxx_fake_c_type, _limb)

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
    template<class Fmpz>
    static typename mp::enable_if<
        traits::is_fmpzxx<Fmpz>, nmodxx_expression>::type red(
                const Fmpz& n, nmodxx_ctx_srcref c)
    {
        return make_nored((n % c.n()).template to<mp_limb_t>(), c);
    }
    template<class Fmpq>
    static typename mp::enable_if<
        traits::is_fmpqxx<Fmpq>, nmodxx_expression>::type red(
                const Fmpq& n, nmodxx_ctx_srcref c)
    {
        return make_nored((n % fmpzxx(c.n())).template to<mp_limb_t>(), c);
    }
    // TODO more

    // only makes sense on immediates
    nmodxx_ctx_srcref _ctx() const {return this->_data().ctx;}
    const nmod_t& _nmod() const {return this->_data().ctx._nmod();}
    void reduce() {NMOD_RED(_limb(), _limb(), _nmod());}
    void set_nored(mp_limb_t n) {this->_data().inner = n;}

    nmodxx_ctx_srcref estimate_ctx() const;

    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx());
    }

    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
};

namespace detail {
struct nmodxx_data;
} // detail

typedef nmodxx_expression<operations::immediate, detail::nmodxx_data> nmodxx;
typedef nmodxx_expression<operations::immediate,
            flint_classes::ref_data<nmodxx, detail::nmodxx_fake_c_type> > nmodxx_ref;
typedef nmodxx_expression<operations::immediate, flint_classes::srcref_data<
    nmodxx, nmodxx_ref, detail::nmodxx_fake_c_type> > nmodxx_srcref;

namespace flint_classes {
template<class Nmod>
struct ref_data<Nmod, detail::nmodxx_fake_c_type>
{
    typedef void IS_REF_OR_CREF;
    typedef Nmod wrapped_t;

    typedef mp_limb_t& data_ref_t;
    typedef const mp_limb_t& data_srcref_t;

    mp_limb_t& inner;
    nmodxx_ctx_srcref ctx;

    ref_data(Nmod& o) : inner(o._data().inner), ctx(o._data().ctx) {}

    static ref_data make(mp_limb_t& f, nmodxx_ctx_srcref ctx)
    {
        return ref_data(f, ctx);
    }

private:
    ref_data(mp_limb_t& fp, nmodxx_ctx_srcref c) : inner(fp), ctx(c) {}
};

template<class Nmod, class Ref>
struct srcref_data<Nmod, Ref, detail::nmodxx_fake_c_type>
{
    typedef void IS_REF_OR_CREF;
    typedef Nmod wrapped_t;

    typedef const mp_limb_t& data_ref_t;
    typedef const mp_limb_t& data_srcref_t;

    const mp_limb_t& inner;
    nmodxx_ctx_srcref ctx;

    srcref_data(const Nmod& o) : inner(o._data().inner), ctx(o._data().ctx) {}
    srcref_data(Ref o) : inner(o._data().inner) {}

    static srcref_data make(const mp_limb_t& f, nmodxx_ctx_srcref ctx)
    {
        return srcref_data(f, ctx);
    }

private:
    srcref_data(const mp_limb_t& fp, nmodxx_ctx_srcref c) : inner(fp), ctx(c) {}
};
} // flint_classes
namespace detail {
struct nmodxx_data
{
    nmodxx_ctx_srcref ctx;
    mp_limb_t inner;
    typedef mp_limb_t& data_ref_t;
    typedef const mp_limb_t& data_srcref_t;

    nmodxx_data(nmodxx_ctx_srcref c) : ctx(c), inner(0) {}

private:
    nmodxx_data(mp_limb_t n, nmodxx_ctx_srcref c)
        : ctx(c), inner(n) {}

public:
    static nmodxx_data make_nored(mp_limb_t n, nmodxx_ctx_srcref c)
    {
        return nmodxx_data(n, c);
    }

    nmodxx_data(const nmodxx_srcref& r)
        : ctx(r.estimate_ctx()), inner(r._limb()) {}
};
} // detail

// Temporary merging isn't really any use here. On the other hand, it does
// not seem to hurt. Let's leave this for now. -- Tom Bachmann (15/10/2013)
#if 0
namespace traits {
template<> struct use_temporary_merging<nmodxx> : mp::false_ { };
} // traits
#endif

namespace traits {
template<class T> struct has_nmodxx_ctx : mp::false_ { };

template<> struct has_nmodxx_ctx<nmodxx> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmodxx_ref> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmodxx_srcref> : mp::true_ { };
} // traits

namespace detail {
struct has_nmodxx_ctx_predicate
{
    template<class T> struct type : traits::has_nmodxx_ctx<T> { };
};

// XXX this is needed for vectors ...
template<class T>
struct get_nmodxx_ctx
{
    static nmodxx_ctx_srcref get(const T& t) {return t._ctx();}
};
template<class T>
nmodxx_ctx_srcref get_nmodxx_ctx_func(const T& t)
{
    return get_nmodxx_ctx<T>::get(t);
}
} // detail
namespace tools {
template<class Expr>
nmodxx_ctx_srcref find_nmodxx_ctx(const Expr& e)
{
    return detail::get_nmodxx_ctx_func(
            tools::find_subexpr<detail::has_nmodxx_ctx_predicate>(e));
}
} // tools
template<class Operation, class Data>
inline nmodxx_ctx_srcref
nmodxx_expression<Operation, Data>::estimate_ctx() const
{
    return tools::find_nmodxx_ctx(*this);
}

namespace traits {
template<class T> struct is_nmodxx : mp::or_<
     traits::is_T_expr<T, nmodxx>,
     flint_classes::is_source<nmodxx, T> > { };
} // traits

namespace rules {
#define NMODXX_COND_S FLINTXX_COND_S(nmodxx)
#define NMODXX_COND_T FLINTXX_COND_T(nmodxx)

#define NMODXX_DEFINE_INSTANTIATE_TEMPORARIES(Classname) \
template<class Expr> \
struct use_default_temporary_instantiation<Expr, Classname> : mp::false_ { }; \
template<class Expr> \
struct instantiate_temporaries<Expr, Classname> \
{ \
    static Classname get(const Expr& e) \
    { \
        return Classname(tools::find_nmodxx_ctx(e)); \
    } \
};

// This is in order to make temporary allocation work even if there is no
// immediate subexpression - c/f test_temporaries
NMODXX_DEFINE_INSTANTIATE_TEMPORARIES(nmodxx)

FLINTXX_DEFINE_EQUALS(nmodxx, e1._limb() == e2._limb())

FLINT_DEFINE_GET_COND(conversion, mp_limb_t, NMODXX_COND_S, from._limb())

template<class Nmod>
struct to_string<Nmod, typename mp::enable_if< NMODXX_COND_S<Nmod> >::type>
{
    static std::string get(const Nmod& i, int base /* ignored */)
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

//////////////////////////////////////////////////////////////////////////////
// NMOD_VEC CLASS AND RULES
//////////////////////////////////////////////////////////////////////////////
namespace detail {
struct nmod_vector_data
{
    slong size;
    mp_limb_t* array;
    nmodxx_ctx_srcref ctx;

    nmod_vector_data(slong n, nmodxx_ctx_srcref c)
        : size(n), array(_nmod_vec_init(n)), ctx(c) {}

    ~nmod_vector_data() {_nmod_vec_clear(array);}

    nmod_vector_data(const nmod_vector_data& o)
        : size(o.size), array(_nmod_vec_init(o.size)), ctx(o.ctx)
    {
        _nmod_vec_set(array, o.array, size);
    }

    nmodxx_ref at(slong i) {return nmodxx_ref::make(array[i], ctx);}
    nmodxx_srcref at(slong i) const {return nmodxx_srcref::make(array[i], ctx);}
};

struct nmod_vector_traits
    : wrapped_vector_traits<nmodxx, slong, nmodxx_ref, nmodxx_srcref, mp_limb_t>
{
    template<class Expr>
    static typename Expr::evaluated_t create_temporary(const Expr& e)
    {
        return typename Expr::evaluated_t(e.size(), tools::find_nmodxx_ctx(e));
    }
};
} // detail

// TODO would it make more sense to have this have its own class?
typedef vector_expression<
    detail::nmod_vector_traits, operations::immediate,
    detail::nmod_vector_data> nmod_vecxx;
// TODO references

namespace traits {
template<> struct has_nmodxx_ctx<nmod_vecxx> : mp::true_ { };
} // traits
namespace detail {
template<>
struct get_nmodxx_ctx<nmod_vecxx>
{
    static nmodxx_ctx_srcref get(const nmod_vecxx& v)
    {
        return v._data().ctx;
    }
};
} // detail

template<>
struct enable_vector_rules<nmod_vecxx> : mp::false_ { };

namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct NMOD_VECXX_COND_S : mp::equal_types<T, nmod_vecxx> { };
#define NMOD_VECXX_COND_T NMOD_VECXX_COND_S

// TODO references
FLINT_DEFINE_GET(equals, bool, nmod_vecxx,
        e1.size() == e2.size()
        && _nmod_vec_equal(e1._data().array, e2._data().array, e1.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, nmod_vecxx,
        NMOD_VECXX_COND_S, NMOD_VECXX_COND_S,
        _nmod_vec_add(to._data().array, e1._data().array, e2._data().array,
            to.size(), to._data().ctx._nmod()))

// TODO more
} // rules
} // flint

#endif
