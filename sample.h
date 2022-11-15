/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_SAMPLE_H
#define _EPP_SAMPLE_H 1

namespace EPP
{
    /**
     * EPP is an algorithim not an implementation, and it can be applied to the clients in memory data
     * therefore it's defined as a template based on the users preferred data model
     **/

    // distributed single precision Measurement value
    class Value
    {
        union
        {
            float floating;
            uint32_t binary;
        };

        inline Value get(uint8_t *&ptr) const noexcept
        {
            Value v;
            v.binary |= *ptr++ << 24; // big endian
            v.binary |= *ptr++ << 16;
            v.binary |= *ptr++ << 8;
            v.binary |= *ptr++;
            return v;
        };

        inline void put(uint8_t *&ptr, Value v) const noexcept
        {
            *ptr++ = v.binary >> 24;
            *ptr++ = v.binary >> 16;
            *ptr++ = v.binary >> 8;
            *ptr++ = v.binary;
        }

        Value(float value) : floating(value){};
        Value(double value) : floating((float)value){};
        Value() : binary(0){};
    };


    template <typename _float>
    class DefaultSample : public Sample
    {
    public:
        inline const double operator()(Event event, Measurement measurement) const noexcept
        {
            return (double)data[Sample::measurements * event + measurement];
        };

        DefaultSample(Measurement measurements,
                      Event events,
                      SampleSubset<DefaultSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        DefaultSample(const Measurement measurements,
                      const Event events,
                      const _float *const data,
                      SampleSubset<DefaultSample> subset) noexcept
            : Sample(measurements, events, subset, NoKey), data(data){};

    protected:
        epp_word get_word(Measurement measurement, Event event) const noexcept
        {
            float f = data[Sample::measurements * event + measurement];
            return *(epp_word *)&f;
        };

        void put_word(Measurement measurement, Event event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[Sample::measurements * event + measurement] = (_float)f;
        };

    private:
        const _float *const data;
    };

    template <typename _float>
    class TransposeSample : public Sample, protected Blob
    {
    public:
        inline double operator()(Event event, Measurement measurement) const noexcept
        {
            return (double)data[Sample::events * measurement + event];
        };

        TransposeSample(Measurement measurements,
                        Event events,
                        SampleSubset<TransposeSample> subset,
                        Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        TransposeSample(const Measurement measurements,
                        const Event events,
                        const _float *const data) noexcept
            : Sample(measurements, events), data(data){};

    protected:
        epp_word get_word(Measurement measurement, Event event) const noexcept
        {
            float f = data[Sample::events * measurement + event];
            return *(epp_word *)&f;
        };

        void put_word(Measurement measurement, Event event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[Sample::events * measurement + event] = (_float)f;
        };

    private:
        const _float *const data;
    };

    template <typename _float>
    class PointerSample : public Sample
    {
    public:
        inline double operator()(Event event, Measurement measurement) const noexcept
        {
            return (double)data[measurement][event];
        };

        PointerSample(Measurement measurements,
                      Event events,
                      SampleSubset<PointerSample> subset,
                      Key key) noexcept
            : Sample(measurements, events, subset, key), data(nullptr){};

        PointerSample(const Measurement measurements,
                      const Event events,
                      const _float *const *const data,
                      SampleSubset<PointerSample> subset) noexcept
            : Sample(measurements, events, subset, NoKey), data(data){};

    protected:
        epp_word get_word(Measurement measurement, Event event) const noexcept
        {
            float f = data[measurement][event];
            return *(epp_word *)&f;
        };

        void put_word(Measurement measurement, Event event, epp_word value) noexcept
        {
            float f = *(float *)&value;
            data[measurement][event] = (_float)f;
        };

    private:
        const _float *const *const data;
    };

    /**
     * communications services assuming control and data streams
     * are already set up. blob faults are handled here.
     * requests and results are safely queued for the pursuer
     **/

    class Remote //: Blob::Handler
    {
        // friend class Blob;

    public:
        enum Service
        {
            request,
            result,
            fault,
            content
        };

        // thread safe send/receive one json message

        void out(const json &encoded); // does not block

        json in(); // blocks calling thread

    protected:
        void copy(
            std::istream *in,
            std::ostream *out,
            std::streamsize count);

        void startFault(
            Key key);

        void transmit();

        void receive();

        Remote();

        ~Remote();

    private:
        std::mutex mutex;
        std::queue<json> incoming;
        std::queue<json> outgoing;
        std::condition_variable wake_in;
        std::condition_variable wake_out;
        std::iostream *remote_control;
        std::iostream *remote_data;
        std::thread receiver;
        std::thread transmitter;
        volatile bool on_the_air = true;
    };
    /**
     * remote worker instance
     **/
    template <class ClientSample>
    class Pursuer;

    typedef DefaultSample<float> CloudSample;

    class CloudPursuer : Remote, public Pursuer<CloudSample>
    {
    public:
        void start(const json &encoded);

        void finish(Request<CloudSample> *request) noexcept;

        void finish(const json &encoded);

        json remote();

        CloudPursuer(
            const Parameters &parameters) noexcept;

        ~CloudPursuer();
    };

    /**
     * utility classes for serializing sample and subset as streams
     */

    class SampleStream : public std::iostream
    {
    protected:
        class sample_buffer : public std::streambuf
        {

        public:
            explicit sample_buffer(Sample &sample);
            virtual ~sample_buffer();

        protected:
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Sample *sample;
            epp_word *buffer;
            long next_event;
        };

    public:
        explicit SampleStream(Sample &sample);
    };

    class SubsetStream : public std::iostream
    {
    protected:
        class subset_buffer : public std::streambuf
        {
        public:
            explicit subset_buffer(Subset &subset);
            virtual ~subset_buffer();
            virtual std::streambuf::int_type underflow();
            virtual std::streambuf::int_type overflow(std::streambuf::int_type value);
            virtual std::streambuf::int_type sync();

        private:
            Subset *subset;
            uint8_t *buffer;
            long next_event;
        };

    public:
        explicit SubsetStream(Subset &subset);
    };
}

#endif /* _EPP_SAMPLE_H */
