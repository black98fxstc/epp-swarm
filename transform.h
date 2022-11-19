
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
    // fftw needs special alignment to take advantage of vector instructions
    class FFTData
    {
        friend class Transform;

        float *data;

    public:
        FFTData()
        {
            data = nullptr;
        }

        ~FFTData()
        {
            if (data)
                fftwf_free(data);
        };

        float *operator*() noexcept
        {
            if (!data)
                data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
            return data;
        };

        inline float &operator[](const int i)
        {
            return data[i];
        }

        void zero() noexcept
        {
            if (!data)
                data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
            std::fill(data, data + (N + 1) * (N + 1), (float)0);
        };

        void dump(const std::string &file)
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
    };

    class Transform
    {
        void *DCT;
        void *IDCT;

    public:
        static void swap(FFTData &red, FFTData &blue)
        {
            float *temp = red.data;
            red.data = blue.data;
            blue.data = temp;
            if (!red.data)
                red.data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
            if (!blue.data)
                blue.data = (float *)fftw_malloc(sizeof(float) * (N + 1) * (N + 1));
        }

        Transform() noexcept
        {
            FFTData in;
            FFTData out;
            // FFTW planning is slow and not thread safe so we do it here
            DCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                            FFTW_REDFT00, FFTW_REDFT00, 0);
            // actually they are the same in this case but leave it for now
            IDCT = (void *)fftwf_plan_r2r_2d((N + 1), (N + 1), *in, *out,
                                             FFTW_REDFT00, FFTW_REDFT00, 0);
            assert(DCT && IDCT);
        };

        ~Transform()
        {
            fftwf_destroy_plan((fftwf_plan)DCT);
            fftwf_destroy_plan((fftwf_plan)IDCT);
        };

        void forward(FFTData &in, FFTData &out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)DCT, *in, *out);
        };

        void reverse(FFTData &in, FFTData &out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)IDCT, *in, *out);
        };
    };
}

#endif /* _TRANSFORM_H */
