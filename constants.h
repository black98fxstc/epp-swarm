#ifndef _EPP_CONSTANTS_H
#define _EPP_CONSTANTS_H 1
namespace EPP
{
    // resolution of the density estimator
    // FFT is fastest when N has lots of small prime factors
    const int N = 1 << 8;
    const double pi = 3.14159265358979323846;

    // incremental width of the kernel relative to full scale, if I've done my sums right
    const double W = .01;
}
#endif /* _EPP_CONSTANTS_H */