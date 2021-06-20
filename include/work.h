#include <client.h>
#include <boundary.h>

#include <queue>

// testing
#include <random>

namespace EPP
{
    // testing stuff
    std::default_random_engine generator;
    std::binomial_distribution<int> binomial(2000, 0.5);
    std::binomial_distribution<int> quality(500, 0.5);
    std::binomial_distribution<int> coin_toss(1, 0.5);

    // abstract class representing a unit of work to be done
    // virtual functions let subclasses specialize tasks
    // handles work_completed and work_outstanding

    std::recursive_mutex mutex;
    std::condition_variable_any work_available;
    std::condition_variable_any work_completed;
    int work_outstanding = 0;

    // these are essential constants that are read only
    // so safely shared by all threads
    struct worker_sample
    {
        const int measurements;
        const long events;
        const float *const data;
        const std::vector<bool> subset;
    };

    class Work
    {
    public:
        const struct worker_sample sample;

        // many threads can execute this in parallel
        virtual void parallel()
        {
            throw std::runtime_error("unimplemented");
        };
        // then only one thread at a time can run
        virtual void serial()
        {
            throw std::runtime_error("unimplemented");
        };

        ~Work()
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            if (--work_outstanding == 0)
                work_completed.notify_all();
        };

    protected:
        explicit Work(const worker_sample &sample) : sample(sample)
        {
            std::unique_lock<std::recursive_mutex> lock(EPP::mutex);
            ++work_outstanding;
        };
    };

    volatile static bool kiss_of_death = false;
    std::queue<Work *> work_list;

    // a generic worker thread. looks for work, does it, deletes it
    // virtual functions in the work object do all the real work
    class Worker
    {
    public:
        Worker()
        {
            std::unique_lock<std::recursive_mutex> lock(mutex);
            while (!kiss_of_death)
                if (work_list.empty())
                    work_available.wait(lock);
                else
                {
                    Work *work = work_list.front();
                    work_list.pop();
                    lock.unlock();
                    work->parallel();
                    lock.lock();
                    work->serial();
                    delete work;
                }
        };
    };

    // pursue a particular X, Y pair
    class PursueProjection : public Work
    {
        // fftw needs special alignment to take advantage of vector instructions
        class FFTData
        {
            float *data;

        public:
            FFTData()
            {
                data = nullptr;
            }

            ~FFTData();

            float *operator*();

            inline float &operator[](const int i)
            {
                return data[i];
            }

            void zero();
        };

        class Transform
        {
            void * DCT;
            void * IDCT;

        public:
            Transform();

            ~Transform();

            void forward(FFTData &in, FFTData &out);

            void reverse(FFTData &in, FFTData &out);
        };

        static Transform transform;

        class Kernel
        {
            float k[N + 1];

        public:
        	const double pi = 3.14159265358979323846;

            void apply(FFTData &cosine, FFTData &filtered, int pass)
            {
            	double width = .001 * N * pass;
				for (int i = 0; i <= N; i++)
					k[i] = exp(- i * i * width * width);

				float *data = *cosine;
				float *smooth = *filtered;
                for (int i = 0; i <= N; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        smooth[i + (N + 1) * j] = data[i + (N + 1) * j] * k[i] * k[j];
						smooth[i + (N + 1) * j] = data[j + (N + 1) * i] * k[j] * k[i];
                    }
					smooth[i + (N + 1) * i] = data[i + (N + 1) * i] * k[i] * k[i];
                }
            }
        };

        static Kernel kernel;

    public:
        const int X, Y;

        std::vector<ColoredEdge<short, bool>> separatrix;
        std::vector<bool> in;
        std::vector<bool> out;

        PursueProjection(
            const worker_sample sample,
            const int X,
            const int Y)
            : Work(sample), X(X), Y(Y){};

        ~PursueProjection()= default;;

        virtual void parallel();

        virtual void serial();
    };

    std::vector<int> qualified_measurements;

    class QualifyMeasurement : public Work
    {
        class Scratch
        {
            float *data;
            long size;

        public:
            Scratch()
            {
                data = nullptr;
                size = 0;
            }

            ~Scratch()
            {
				delete[] data;
            }

            float *&reserve(long size)
            {
                if (this->size < size)
                {
                    delete[] data;
                    data = nullptr;
                }
                if (!data)
                {
                    data = new float[size];
                    this->size = size;
                }
                return data;
            }

            inline float &operator[](const int i)
            {
                return data[i];
            }
        };

        static thread_local Scratch scratch;

    public:
        const int X;
        double KLDn = 0;
        double KLDe = 0;
        bool qualified = false;

        QualifyMeasurement(
            const worker_sample sample,
            const int X)
            : Work(sample), X(X){};

        virtual void parallel();

        virtual void serial();
    };
}