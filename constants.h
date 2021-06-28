namespace EPP
{
    // resolution of the density estimator
    // FFT is fastest when N has lots of small prime factors
    const int N = 1 << 8;
    // incremental width of the kernel relative to full scale // if I've done my sums right
    const double W = .01;
}