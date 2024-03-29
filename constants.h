
/*
 * Developer: Wayne Moore <wmoore@stanford.edu> 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_CONSTANTS_H
#define _EPP_CONSTANTS_H 1

#include <cstdint>

namespace EPP
{
    // resolution of the density estimator
    // FFT is fastest when N has lots of small prime factors
    const int N = 1 << 8;

    const unsigned int max_passes = 10;

    const double pi = 3.14159265358979323846;
    const double sqrt2 = sqrt(2);
    const double sqrt_pi = sqrt(pi);

    // random data for the non-key value
    const std::uint8_t no_key[32] =
        {
            248, 150, 225, 10, 94, 55, 71, 31,
            181, 210, 149, 28, 122, 48, 112, 243,
            61, 249, 88, 160, 217, 176, 40, 120,
            186, 219, 232, 211, 216, 152, 239, 141};

    extern int verbose_flags;
}
#endif /* _EPP_CONSTANTS_H */
