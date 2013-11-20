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

#ifndef FMPZ_FACTORXX_H
#define FMPZ_FACTORXX_H

#include "fmpz_factor.h"
#include "fmpz_vec.h"

#include "flintxx/ltuple.h"

// TODO codegen
// TODO factor_pp1 multiple return values

namespace flint {
FLINT_DEFINE_THREEARY(factor_trial_range)
FLINT_DEFINE_UNOP(expand)
FLINT_DEFINE_UNOP(expand_iterative)
FLINT_DEFINE_UNOP(expand_multiexp)

namespace detail {
template<class Delay>
class fmpz_factorxx_delayed
{
private:
    fmpz_factor_t inner;

    void copy_init(const fmpz_factorxx_delayed& o)
    {
	_fmpz_factor_fit_length(inner, o.inner->num);
	_fmpz_factor_set_length(inner, o.inner->num);
	inner->sign = o.inner->sign;
	for(slong i = 0;i < o.inner->num;++i)
	{
	    fmpz_set(inner->p + i, o.inner->p + i);
	    inner->exp[i] = o.inner->exp[i];
	}
    }

public:
    fmpz_factorxx_delayed() {fmpz_factor_init(inner);}
    ~fmpz_factorxx_delayed() {fmpz_factor_clear(inner);}

    fmpz_factorxx_delayed(const fmpz_factorxx_delayed& o)
    {
	fmpz_factor_init(inner);
	copy_init(o);
    }

    fmpz_factorxx_delayed& operator=(const fmpz_factorxx_delayed& o)
    {
	copy_init(o);
	return *this;
    }

    bool operator==(const fmpz_factorxx_delayed& o)
    {
	if(o.sign() != sign() || o.size() != size())
	    return false;
	for(ulong i = 0;i < size();++i)
	    if(p(i) != o.p(i) || exp(i) != o.exp(i))
		return false;
	return true;
    }

    ulong size() const {return inner->num;}
    ulong exp(slong i) const {return inner->exp[i];}
    ulong& exp(slong i) {return inner->exp[i];}
    fmpzxx_srcref p(slong i) const {return fmpzxx_srcref::make(inner->p + i);}
    fmpzxx_ref p(slong i) {return fmpzxx_ref::make(inner->p + i);}
    int sign() const {return inner->sign;}
    int& sign() {return inner->sign;}

    fmpz_factor_t& _data() {return inner;}
    const fmpz_factor_t& _data() const {return inner;}

    void print() const {fmpz_factor_print(inner);}

    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_factor(const Fmpz& f)
    {
	fmpz_factor(_data(), f.evaluate()._fmpz());
    }

    template<class T>
    typename mp::enable_if<traits::fits_into_slong<T> >::type
    set_factor(T t)
    {
	fmpz_factor_si(_data(), t);
    }

    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz>, bool>::type
    set_factor_trial_range(const Fmpz& f, ulong start, ulong nprimes)
    {
	return fmpz_factor_trial_range(_data(), f.evaluate()._fmpz(),
		start, nprimes);
    }

    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz>, bool>::type
    set_factor_pp1(const Fmpz& f, ulong B1, ulong B2_sqrt, ulong c)
    {
	return fmpz_factor_pp1(_data(), f.evaluate()._fmpz(),
		B1, B2_sqrt, c);
    }

#define FLINTXX_DEFINE_MEMBER_UNOP_EXTRA(funcname, Class, rtype) \
    FLINT_UNOP_BUILD_RETTYPE(funcname, rtype, Class) \
    funcname() const {return flint::funcname(*this);}

    FLINTXX_DEFINE_MEMBER_UNOP_EXTRA(expand, fmpz_factorxx_delayed, fmpzxx)
    FLINTXX_DEFINE_MEMBER_UNOP_EXTRA(expand_iterative, fmpz_factorxx_delayed, fmpzxx)
    FLINTXX_DEFINE_MEMBER_UNOP_EXTRA(expand_multiexp, fmpz_factorxx_delayed, fmpzxx)
};
} // detail
typedef detail::fmpz_factorxx_delayed<void> fmpz_factorxx;

template<class Fmpz>
inline typename mp::enable_if<mp::or_<
        traits::is_fmpzxx<Fmpz>, traits::fits_into_slong<Fmpz> >,
    fmpz_factorxx>::type factor(const Fmpz& f)
{
    fmpz_factorxx res;
    res.set_factor(f);
    return res;
}

namespace rules {
namespace rdetail {
typedef make_ltuple<mp::make_tuple<bool, fmpz_factorxx>::type>::type
    fmpz_factor_rt;

template<class T> struct signed_or_fmpz
    : mp::or_<FMPZXX_COND_S<T>, traits::fits_into_slong<T> > { };
} // rdetail
FLINT_DEFINE_THREEARY_EXPR_COND3(factor_trial_range_op, rdetail::fmpz_factor_rt,
	rdetail::signed_or_fmpz,
	traits::is_unsigned_integer, traits::is_unsigned_integer,
	to.template get<0>() = to.template get<1>().set_factor_trial_range(
	    e1, e2, e3))

FLINT_DEFINE_UNARY_EXPR_(expand_op, fmpzxx, fmpz_factorxx,
	fmpz_factor_expand(to._fmpz(), from._data()))
FLINT_DEFINE_UNARY_EXPR_(expand_iterative_op, fmpzxx, fmpz_factorxx,
	fmpz_factor_expand_iterative(to._fmpz(), from._data()))
FLINT_DEFINE_UNARY_EXPR_(expand_multiexp_op, fmpzxx, fmpz_factorxx,
	fmpz_factor_expand_multiexp(to._fmpz(), from._data()))
} // rules

template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, fmpz_factorxx>::type
factor_pp1(const Fmpz& f, ulong B1, ulong B2_sqrt, ulong c)
{
    fmpz_factorxx res;
    res.set_factor_pp1(f, B1, B2_sqrt, c);
    return res;
}

inline void print(const fmpz_factorxx& f)
{
    f.print();
}
} // flint

#endif
