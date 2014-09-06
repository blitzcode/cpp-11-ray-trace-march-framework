
#ifndef TRACE_H
#define TRACE_H

#include <string>

#include "types.h"

void Trace(const char *fmt, ...);
std::string PrintBytesHumanReadable(uint64 bytes);

#endif // TRACE_H

