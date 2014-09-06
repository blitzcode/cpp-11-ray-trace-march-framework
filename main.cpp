
#include "application.h"
#include "timer.h"
#include "sampling.h"

int main(int argc, char **argv)
{
    // Initialize timer while we're still single threaded
    TimerGetTick();

    SAMP::Initialize();

    Application app;
    return app.Main(argc, argv);
}

