
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
        const unsigned N;
        void *DCT;
        void *IDCT;
        std::recursive_mutex &mutex;

    public:
        class Data
        {
            Transform &transform;
        
        public:
            float *data;

            inline void clear ()
            {
                std::fill(data, data + (transform.N + 2) * (transform.N + 2), (float)0);
            }

            inline float &operator[] (const int i)
            {
                return data[i];
            }

            Data (Transform &transform) : transform(transform)
            {
                std::lock_guard<std::recursive_mutex> lock(transform.mutex);
                data = (float *)fftw_malloc(sizeof(float) * (transform.N + 2) * (transform.N + 2));
            }

            ~Data ()
            {
                std::lock_guard<std::recursive_mutex> lock(transform.mutex);
                fftwf_free(data);
            }
        };

        // FFTW not thread safe except execute
        // this takes a long time so do it in the constructor
        Transform(unsigned N, std::recursive_mutex &mutex) : N(N), mutex(mutex)
        {
            std::lock_guard<std::recursive_mutex> lock(mutex);
            float *in = (float *)fftw_malloc(sizeof(float) * (N + 2) * (N + 2));
            float *out = (float *)fftw_malloc(sizeof(float) * (N + 2) * (N + 2));
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
            std::lock_guard<std::recursive_mutex> lock(mutex);
            fftwf_destroy_plan((fftwf_plan)DCT);
            fftwf_destroy_plan((fftwf_plan)IDCT);
        };

        // thread safe read only access to plan, pointers reference thread_local data
        void forward(Data &in, Data &out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)DCT, in.data, out.data);
        };

        void reverse(Data &in, Data &out) noexcept
        {
            fftwf_execute_r2r((fftwf_plan)IDCT, in.data, out.data);
        };

        void dump(Data &data, const std::string &file)
        {
            std::ofstream out(file, std::ios::out);
            for (unsigned i = 0; i < (N + 1) * (N + 1);)
            {
                out << data.data[i++];
                for (int j = N; j > 0; --j)
                    out << "," << data[i++];
                out << std::endl;
            }
            out.close();
        }
    };
}

#endif /* _TRANSFORM_H */
