
#include "timer.h"

#include <chrono>

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

