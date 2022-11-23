
/*
 * Developer: Wayne Moore <wmoore@stanford.edu> 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _TRANSFORM_H
#define _TRANSFORM_H 1

#include <iostream>

#include <fftw3.h>

namespace EPP
{

    class Transform
    {
        const unsigned  N;
        void *DCT;
        void *IDCT;

    public:

        Transform(unsigned N, std::recursive_mutex &mutex) : N(N), mutex(mutex)
        {
            float* in = allocate();
            float* out = allocate();
            // FFTW planning is slow and not thread safe so we do it here
            DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), in, out,
                                            FFTW_REDFT00, FFTW_REDFT00, 0);
            // actually they are the same in this case but leave it for now
            IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), in, out,
                                             FFTW_REDFT00, FFTW_REDFT00, 0);
            assert(DCT && IDCT);
            fftwf_free(in);
            fftwf_free(out);
        };

        ~Transform()
        {
            fftwf_destroy_plan((fftwf_plan)DCT);
            fftwf_destroy_plan((fftwf_plan)IDCT);
            for (float *fft_data : allocated)
                fftwf_free(fft_data);
        };

        void allocate(float* &fft_data)
        {
            if (fft_data)
                return;
            fft_data = allocate();
            std::lock_guard<std::recursive_mutex> lock(mutex);
            this->allocated.push_back(fft_data);
        }

        void forward(float *in, float *out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)DCT, in, out);
        };

        void reverse(float *in, float *out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)IDCT, in, out);
        };

        void dump(float *data, const std::string &file)
        {
            std::ofstream out(file, std::ios::out);
            for (int i = 0; i <(N + 1) * (N + 1); )
            {
                out << data[i++];
                for (int j = N; j > 0; --j)
                    out << "," << data[i++];
                out << std::endl;
            }
            out.close();
        }

    protected:
        std::vector<float *> allocated;
        std::recursive_mutex &mutex;

        float* allocate()
        {
            return (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        }
    };
}

#endif /* _TRANSFORM_H */
