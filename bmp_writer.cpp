
#include "bmp_writer.h"

#include <cstdio>
#include <cassert>

struct __attribute__ ((packed)) BmpHeader
{
    uint16 type;        // Must be 'BM'
    uint32 size_file;
    uint16 reserved1;
    uint16 reserved2;
    uint32 offs_bits;   // In bytes
    uint32 bmih_size;
    int32  width;
    int32  height;
    uint16 planes;      // Must be 1
    uint16 bitcount;
    uint32 compression; // BI_RGB = 0
    uint32 size_image;  // All below set to 0
    int32  xpelsperm;
    int32  ypelsperm;
    uint32 clr_used;
    uint32 clr_imp;
};

void WriteBitmap(const char *filename, uint width, uint height, uint32 *bitmap)
{
    BmpHeader header;
    header.type        = (uint16('M') << 8) + uint16('B');
    header.size_file   = sizeof(BmpHeader) + (width * height * 4);
    header.reserved1   = 0;
    header.reserved2   = 0;
    header.offs_bits   = sizeof(BmpHeader);
    header.bmih_size   = header.offs_bits - 14;
    header.width       = width;
    header.height      = height;
    header.planes      = 1;
    header.bitcount    = 32;
    header.compression = 0;
    header.size_image  = 0;
    header.xpelsperm   = 0;
    header.ypelsperm   = 0;
    header.clr_used    = 0;
    header.clr_imp     = 0;

    std::FILE *file = std::fopen(filename, "wb");
    assert(file != nullptr);
    size_t ret;
    ret = std::fwrite(&header, sizeof(BmpHeader), 1, file);
    assert(ret == 1);
    ret = std::fwrite(bitmap, 4, width * height, file);
    assert(ret == width * height);
    ret = std::fclose(file);
    assert(ret == 0);
}

