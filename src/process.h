#pragma once

#ifdef WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace tmd {

inline const int my_pid() {
#ifdef WIN32
    return GetCurrentProcessId();
#else
    return getpid();
#endif
}

}