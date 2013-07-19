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

#ifndef FLINTXX_FLINT_EXCEPTION_H
#define FLINTXX_FLINT_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace flint {
class flint_exception
    : public std::domain_error // ?
{
public:
    flint_exception(const std::string& what)
        : std::domain_error("FLINT: " + what) {}
};

void execution_check(bool worked, const std::string& where,
        const std::string& context)
{
    if(!worked)
        throw flint_exception(context + " computation failed: " + where);
}
}

#endif
