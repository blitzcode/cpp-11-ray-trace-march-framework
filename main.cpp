
#include "application.h"
#include "timer.h"

int main(int argc, char **argv)
{
    // Initialize timer while we're still single threaded
    TimerGetTick();

    Application app;
    return app.Main(argc, argv);
}

