
#include "timer.h"

#include <chrono>
#include <ctime>

std::string DateTimeString()
{
    // Return a filename friendly date + time string

    auto now       = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm *tm    = std::localtime(&in_time_t);

    char buf[128];
    std::snprintf(
        buf,
        sizeof(buf),
        "%i-%i-%i_%ih-%im-%is",
        1900 + tm->tm_year,
        tm->tm_mday,
        tm->tm_mon + 1,
        tm->tm_hour,
        tm->tm_min,
        tm->tm_sec);

    return std::string(buf);
}

double TimerGetTick()
{
    // High precision timer, return in seconds

    // Make returned tick smaller
    static bool first_call = true;
    static std::chrono::time_point<std::chrono::high_resolution_clock> start;
    if (first_call)
    {
        first_call = false;
        start = std::chrono::high_resolution_clock::now();
    }

    return double(std::chrono::duration_cast<std::chrono::microseconds>
        (std::chrono::high_resolution_clock::now() - start).count())
        / 1000.0 / 1000.0;
}

