# Try to find the PThreads librairies
# PThreads_FOUND - system has PThreads lib
# PThreads_INCLUDE_DIRS - the PThreads include directory
# PThreads_LIBRARIES - Libraries needed to use PThreads

if (PThreads_INCLUDE_DIRS AND PThreads_LIBRARIES)
		# Already in cache, be silent
		set(PThreads_FIND_QUIETLY TRUE)
endif (PThreads_INCLUDE_DIRS AND PThreads_LIBRARIES)

find_path(PThreads_INCLUDE_DIRS NAMES pthread.h )
find_library(PThreads_LIBRARIES NAMES pthreads libpthreads )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PThreads DEFAULT_MSG PThreads_INCLUDE_DIRS PThreads_LIBRARIES)

mark_as_advanced(PThreads_INCLUDE_DIRS PThreads_LIBRARIES)