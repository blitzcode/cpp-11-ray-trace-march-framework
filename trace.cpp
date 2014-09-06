
#include "trace.h"

#include <cstdio>
#include <cstdarg>
#include <thread>
#include <string>
#include <sstream>

#include "timer.h"

void Trace(const char *fmt, ...)
{
    // Lead with tick and thread ID
    const double tick = TimerGetTick();
    std::stringstream stream;
    stream << std::this_thread::get_id();
    std::printf("%-14s @ %6.2fs - ", stream.str().c_str(), tick);

    // Passed format string
	va_list argp;
	va_start(argp, fmt);
    char fmt_buf[1024];
    std::vsnprintf(fmt_buf, sizeof(fmt_buf), fmt, argp);
	va_end(argp);
    std::printf("%s\n", fmt_buf);
}

