
#include "types.h"

void TestTypes()
{
    // Static size checks
    static_assert(sizeof(int8)   == 1, "int8 size test failed"  );
    static_assert(sizeof(int16)  == 2, "int16 size test failed" );
    static_assert(sizeof(int32)  == 4, "int32 size test failed" );
    static_assert(sizeof(int64)  == 8, "int64 size test failed" );
    static_assert(sizeof(uint8)  == 1, "uint8 size test failed" );
    static_assert(sizeof(uint16) == 2, "uint16 size test failed");
    static_assert(sizeof(uint32) == 4, "uint32 size test failed");
    static_assert(sizeof(uint64) == 8, "uint64 size test failed");
}

