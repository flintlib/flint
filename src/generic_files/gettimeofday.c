/* We only want this file if someone is using pure windows, and not MinGW */

#if !defined(__MINGW64__) && !defined(__MINGW32__) && !defined(__CYGWIN__) && !defined(__CYGWIN32__) && (defined(_WIN64) || defined(_WIN32) || defined(_MSC_VER))

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <time.h>

#include "gettimeofday.h"

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;
    static int      tzflag;

    if(tv)
    {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;
        t -= EPOCHFILETIME;
        t /= 10;
        tv->tv_sec  = (long)(t / 1000000);
        tv->tv_usec = (long)(t % 1000000);
    }

    if (tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;
}

#else

typedef int this_file_is_empty;

#endif /* defined(_WIN64) */
