/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINTXX_FLINT_EXCEPTION_H
#define FLINTXX_FLINT_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace flint {
// This is the common flintxx exception class, which is raised whenever the
// C++ layer diagnoses a problem.
// Note that the C layer will sometimes flint_abort() with a message (in this case
// no exception is thrown).
class flint_exception
    : public std::domain_error // ?
{
public:
    flint_exception(const std::string& what)
        : std::domain_error("FLINT: " + what) {}
};

// Helper function. If worked is true, does nothing, else raises
// flint_exception.
inline void execution_check(bool worked, const std::string& where,
        const std::string& context)
{
    if(!worked)
        throw flint_exception(context + " computation failed: " + where);
}
} // flint

#endif
