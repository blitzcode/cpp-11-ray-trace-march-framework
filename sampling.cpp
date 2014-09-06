
#include "sampling.h"
#include <algorithm>
#include <cmath>

// Some references
//
// - PBRT Chapter 7, Sampling & Reconstruction
//   http://graphics.stanford.edu/~mmp/chapters/pbrt_chapter7.pdf
//
// - Good permutations for scrambled Halton sequences in terms of L2-discrepancy
//   http://www.cs.kuleuven.ac.be/publicaties/rapporten/tw/TW406.pdf
//
// - Random and Deterministic Digit Permutations of the Halton Sequence
//   http://www.math.fsu.edu/~goncharo/Papers/RandomDigitPerm%20v2.pdf
//
// - Efficient Multidimensional Sampling
//   http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.14.9428&rep=rep1&type=pdf

namespace SAMP
{

void InitBWTable();
void InitFaureTable();
void InitReverseTable();
void InitRandomizedTable();

void Initialize()
{
    // Initialize the sampling module
    InitBWTable();
    InitFaureTable();
    InitReverseTable();
    InitRandomizedTable();
}

uint PrimeToIndex(uint prime)
{
    // Get the index of a prime number

    const uint *idx = std::lower_bound(g_prime_table, g_prime_table + PRIME_TBL_SIZE, prime);

    if ((* idx) == prime)
        return idx - g_prime_table;
    else
        return std::numeric_limits<uint>::max();
}

// TODO: Would probably be best to convert these to linear arrays + index array
std::vector<uint> g_braaten_weller_table[BW_TBL_SIZE];
std::vector<uint> g_faure_table[FAURE_TBL_SIZE];
std::vector<uint> g_reverse_table[REVERSE_TBL_SIZE];
std::vector<uint> g_randomized_table[RANDOMIZED_TBL_SIZE];

void InitBWTable()
{
    static_assert(BW_TBL_SIZE <= PRIME_TBL_SIZE, "table size error");

    // Initialize the Braaten-Weller permutation table for the first BW_TBL_SIZE
    // primes

    // This table is from the technical report of Vandewoestyne / Cools, not the
    // original paper. There seems to be no larger table available
    const uint bw_permut_16[] =
    {
        0, 1, 0, 2, 1, 0, 2, 4, 1, 3, 0, 3, 5, 1, 6, 2, 4, 0, 5, 8, 2, 10, 3, 6, 1, 9, 4,
        7, 0, 6, 10, 2, 8, 4, 12, 1, 9, 5, 11, 3, 7, 0, 8, 13, 3, 11, 5, 16, 1, 10, 7, 14,
        4, 12, 2, 15, 6, 9, 0, 9, 14, 3, 17, 6, 11, 1, 15, 7, 12, 4, 18, 8, 2, 16, 10, 5,
        13, 0, 11, 17, 4, 20, 7, 13, 2, 22, 9, 15, 5, 18, 1, 14, 10, 21, 6, 16, 3, 19, 8,
        12, 0, 14, 22, 5, 18, 9, 27, 2, 20, 11, 25, 7, 16, 3, 24, 13, 19, 6, 28, 10, 1, 23,
        15, 12, 26, 4, 17, 8, 21, 0, 16, 8, 26, 4, 22, 13, 29, 2, 19, 11, 24, 6, 20, 14,
        28, 1, 17, 9, 30, 10, 23, 5, 21, 15, 3, 27, 12, 25, 7, 18, 0, 18, 28, 6, 23, 11,
        34, 3, 25, 14, 31, 8, 20, 36, 1, 16, 27, 10, 22, 13, 32, 4, 29, 17, 7, 35, 19, 2,
        26, 12, 30, 9, 24, 15, 33, 5, 21, 0, 20, 31, 7, 26, 12, 38, 3, 23, 34, 14, 17, 29,
        5, 40, 10, 24, 1, 35, 18, 28, 9, 33, 15, 21, 4, 37, 13, 30, 8, 39, 19, 25, 2, 32,
        11, 22, 36, 6, 27, 16, 0, 21, 32, 7, 38, 13, 25, 3, 35, 17, 28, 10, 41, 5, 23, 30,
        15, 37, 1, 19, 33, 11, 26, 42, 8, 18, 29, 4, 39, 14, 22, 34, 6, 24, 12, 40, 2, 31,
        20, 16, 36, 9, 27, 0, 23, 35, 8, 41, 14, 27, 3, 44, 18, 31, 11, 37, 5, 25, 39, 16,
        21, 33, 1, 46, 12, 29, 19, 42, 7, 28, 10, 36, 22, 4, 43, 17, 32, 13, 38, 2, 26, 45,
        15, 30, 6, 34, 20, 40, 9, 24, 0, 26, 40, 9, 33, 16, 49, 4, 36, 21, 45, 12, 29, 6,
        51, 23, 38, 14, 43, 1, 30, 19, 47, 10, 34, 24, 42, 3, 27, 52, 15, 18, 39, 7, 46,
        22, 32, 5, 48, 13, 35, 25, 8, 44, 31, 17, 50, 2, 37, 20, 28, 11, 41
    };

    uint offs = 0;
    for (uint i=0; i<BW_TBL_SIZE; i++)
    {
        const uint base = g_prime_table[i];
        g_braaten_weller_table[i].resize(base);
        for (uint digit=0; digit<base; digit++)
        {
            assert(bw_permut_16[offs + digit] < base);
            g_braaten_weller_table[i][digit] = bw_permut_16[offs + digit];
        }
        offs += base;
    }
    assert(offs == sizeof(bw_permut_16) / sizeof(uint));
}

void GenerateFaurePermutations(uint base, uint *perm_out)
{
    // See Keller's 'Monte Carlo and Beyond' slide 96 or
    // 'Production Rendering' p225 by Ian Stephenson (Ed.)

    assert(base >= 2);

    if (base == 2) // Identity
    {
        perm_out[0] = 0;
        perm_out[1] = 1;
    }
    else if (base % 2 == 0) // Even
    {
        GenerateFaurePermutations(base / 2, perm_out);
        GenerateFaurePermutations(base / 2, &perm_out[base / 2]);
        for (uint i=0; i<base / 2; i++)
            perm_out[i] = 2 * perm_out[i];
        for (uint i=base / 2; i<base; i++)
            perm_out[i] = 2 * perm_out[i] + 1;
    }
    else // Odd
    {
        GenerateFaurePermutations(base - 1, perm_out);
        const uint bm1d2 = (base - 1) / 2;
        for (uint i=0; i<base - 1; i++)
            if (perm_out[i] >= bm1d2)
                perm_out[i]++;
        for (uint i=base - 1; i>=bm1d2 + 1; --i)
            perm_out[i] = perm_out[i - 1];
        perm_out[bm1d2] = bm1d2;
    }
}

void InitFaureTable()
{
    static_assert(FAURE_TBL_SIZE <= PRIME_TBL_SIZE, "table size error");

    // Initialize the Faure permutation table for the first FAURE_TBL_SIZE
    // primes
    for (uint i=0; i<FAURE_TBL_SIZE; i++)
    {
        g_faure_table[i].resize(g_prime_table[i]);
        GenerateFaurePermutations(g_prime_table[i], &g_faure_table[i][0]);
    }
}

void InitReverseTable()
{
    // Reverse permutation by Vandewoestyne & Cools. This does actually work
    // nowhere near as good as the original paper suggests, later papers agree
    // with that assessment. Barely better than the unscrambled Halton
    // sequence

    static_assert(REVERSE_TBL_SIZE <= PRIME_TBL_SIZE, "table size error");

    // It seems almost pointless to store this in a table, but it simplifies
    // the code
    for (uint i=0; i<REVERSE_TBL_SIZE; i++)
    {
        const uint base = g_prime_table[i];
        g_reverse_table[i].resize(base);
        for (uint digit=0; digit<base; digit++)
        {
            if (digit == 0)
                g_reverse_table[i][0] = 0;
            else
                g_reverse_table[i][digit] = base - digit;
        }
    }
}

void InitRandomizedTable()
{
    // Randomized digit permutations

    static_assert(RANDOMIZED_TBL_SIZE <= PRIME_TBL_SIZE, "table size error");

    for (uint i=0; i<RANDOMIZED_TBL_SIZE; i++)
    {
        const uint base = g_prime_table[i];
        g_randomized_table[i].resize(base);
        for (uint digit=0; digit<base; digit++)
            g_randomized_table[i][digit] = digit;

        if (base == 2)
            continue;

        // Knuth Shuffle (http://en.wikipedia.org/wiki/Knuth_shuffle)
        // TODO: Use Mersenne Twister as RNG here
        std::random_shuffle(&g_randomized_table[i][0], &g_randomized_table[i][0] + base);
    }
}

double RadicalInverse(uint n, uint base, const uint *perm)
{
    const double inv_base = 1.0 / double(base);
    double inv_base_i = inv_base;
    double val = 0.0;

    while (n > 0)
    {
        unsigned int digit = (n % base);
        digit = (perm == nullptr) ? digit : perm[digit];
        val += digit * inv_base_i;
        inv_base_i *= inv_base;
        n /= base;
    }

    return val;
}

double RadicalInverseBase2(uint32 n, uint32 scramble)
{
    // Optimized implementation for base 2

    // This just reverses the digits by swapping them with
    // smaller and smaller intervals
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);

    // Random digit scrambling in base 2 is just a XOR
    n ^= scramble;

    // The +1 is to stay in our [0, 1) interval
    return double(n) / double(0xFFFFFFFFULL + 1ULL);
}

double SobolRadicalInverseBase2(uint32 n, uint32 scramble)
{
    for (uint32 v=1 << 31; n!=0; n>>=1, v^=v >> 1)
    {
        if (n & 0x1)
            scramble ^= v;
    }
    return double(scramble) / double(0xFFFFFFFFULL + 1ULL);
}

double LarcherPillichshammerRadicalInverseBase2(uint32 n, uint32 scramble)
{
    for (uint32 v=1 << 31; n; n>>=1, v|=v >> 1)
    {
        if (n & 0x1)
            scramble ^= v;
    }
    return double(scramble) / double(0xFFFFFFFFULL + 1ULL);
}

double FoldedRadicalInverse(uint n, uint base)
{
    const double inv_base = 1.0 / double(base);
    double inv_base_i = inv_base;
    double val = 0.0;
    uint offset = 0;

    while (n + base * inv_base_i != n)
    {
        uint digit = (n + offset) % base;
        val += digit * inv_base_i;
        inv_base_i *= inv_base;
        n /= base;
        offset++;
    }

    return val;
}

double HaltonZarembaSequence(uint n, uint dim)
{
    return FoldedRadicalInverse(n, g_prime_table[dim]);
}

double HammersleyZarembaSequence(uint n, uint dim, uint num_smp)
{
    if (dim == 0)
        return float(n) / float(num_smp);
    else
        return FoldedRadicalInverse(n, g_prime_table[dim - 1]);
}

double CranleyPattersonRotation(double x, double e)
{
    double ret = x + e;
    if (x > 1.0)
        return ret - 1.0;
    else
        return ret;
}

} // namespace SAMP

