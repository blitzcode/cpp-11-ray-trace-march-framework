
#include "trace.h"

#include <cstdio>
#include <cstdarg>
#include <thread>
#include <sstream>

#include "timer.h"

void Trace(const char *fmt, ...)
{
    // Lead with tick and thread ID
    const double tick = TimerGetTick();
    std::stringstream stream;
    stream << std::this_thread::get_id();

    // Passed format string
	va_list argp;
	va_start(argp, fmt);
    char fmt_buf[2048];
    std::vsnprintf(fmt_buf, sizeof(fmt_buf), fmt, argp);
	va_end(argp);

    std::printf("%-14s @ %6.2fs - %s\n", stream.str().c_str(), tick, fmt_buf);
}

std::string PrintBytesHumanReadable(uint64 bytes)
{
    // Convert byte count into more compact and readable KB / MB / GB string

    char buf[64];

    if (bytes < 1024)
        std::snprintf(buf, sizeof(buf), "%iB", int(bytes));
    else if (bytes < 1024 * 1024)
        std::snprintf(buf, sizeof(buf), "%.1fKB", float(double(bytes) / 1024.0));
    else if (bytes < 1024 * 1024 * 1024)
        std::snprintf(buf, sizeof(buf), "%.1fMB", float(double(bytes) / 1024.0 / 1024.0));
    else
        std::snprintf(buf, sizeof(buf), "%.2fGB", float(double(bytes) / 1024.0 / 1024.0 / 1024.0));

    return std::string(buf);
}

