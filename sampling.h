
#ifndef SAMPLING_H
#define SAMPLING_H

#include <cstddef>
#include <vector>
#include <cassert>

#include "types.h"

namespace SAMP
{

// Global initialization function for the sampling module
void Initialize();

// Prime numbers
const size_t PRIME_TBL_SIZE = 1000;
extern const uint g_prime_table[PRIME_TBL_SIZE];
uint PrimeToIndex(uint prime);

// Radical inverse functions
double RadicalInverse(uint n, uint base, const uint *perm = NULL);
double FoldedRadicalInverse(uint n, uint base);
double RadicalInverseBase2(uint32 n, uint32 scramble = 0);
double SobolRadicalInverseBase2(uint32 n, uint32 scramble = 0);
double LarcherPillichshammerRadicalInverseBase2(uint32 n, uint32 scramble = 0);

// Halton & Hammersley sequences
template <class S> double HaltonSequence(uint n, uint dim);
template <class S> double HammersleySequence(uint n, uint dim, uint num_smp);
double HaltonZarembaSequence(uint n, uint dim);
double HammersleyZarembaSequence(uint n, uint dim, uint num_smp);

double CranleyPattersonRotation(double x, double e);

// Braaten-Weller permutations
const size_t BW_TBL_SIZE = 16;
extern std::vector<uint> g_braaten_weller_table[BW_TBL_SIZE];
struct ScrambleBraatenWeller
{
    const uint * operator () (uint prime_base_idx) const
    {
        assert(g_braaten_weller_table[0].empty() == false);
        if (prime_base_idx < BW_TBL_SIZE)
            return &g_braaten_weller_table[prime_base_idx][0];
        else
            return NULL;
    }
};

// Faure permutations
const size_t FAURE_TBL_SIZE = 128;
extern std::vector<uint> g_faure_table[FAURE_TBL_SIZE];
struct ScrambleFaure
{
    const uint * operator () (uint prime_base_idx) const
    {
        assert(g_faure_table[0].empty() == false);
        if (prime_base_idx < FAURE_TBL_SIZE)
            return &g_faure_table[prime_base_idx][0];
        else
            return NULL;
    }
};

// Reverse permutations
const size_t REVERSE_TBL_SIZE = 128;
extern std::vector<uint> g_reverse_table[REVERSE_TBL_SIZE];
struct ScrambleReverse
{
    const uint * operator () (uint prime_base_idx) const
    {
        assert(g_reverse_table[0].empty() == false);
        if (prime_base_idx < REVERSE_TBL_SIZE)
            return &g_reverse_table[prime_base_idx][0];
        else
            return NULL;
    }
};

// Randomized permutations
const size_t RANDOMIZED_TBL_SIZE = 128;
extern std::vector<uint> g_randomized_table[RANDOMIZED_TBL_SIZE];
struct ScrambleRandomized
{
    const uint * operator () (uint prime_base_idx) const
    {
        assert(g_randomized_table[0].empty() == false);
        if (prime_base_idx < RANDOMIZED_TBL_SIZE)
            return &g_randomized_table[prime_base_idx][0];
        else
            return NULL;
    }
};

// No scrambling
struct ScrambleNone
{
    const uint * operator () (uint) const { return NULL; }
};

//
// Template inline implementations
//

template <class S> double HaltonSequence(uint n, uint dim)
{
    S scramble;
    return RadicalInverse(n, g_prime_table[dim], scramble(dim));
}

template <class S> double HammersleySequence(uint n, uint dim, uint num_smp)
{
    S scramble;
    if (dim == 0)
        return double(n) / double(num_smp);
    else
        return RadicalInverse(n, g_prime_table[dim - 1], scramble(dim));
}

} // namespace SAMP

#endif // SAMPLING_H

