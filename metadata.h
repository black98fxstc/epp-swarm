
/*
 * Developer: Wayne Moore <wmoore@stanford.edu>
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_METADATA_H
#define _EPP_METADATA_H 1

/**
 * Provides a content based associative memory service
 * of blobs shared between clients and servers
 **/

struct Meta
{
    std::iostream *stream = nullptr;
    std::uint8_t *_buffer = nullptr;
    unsigned long int size = 0;
    bool valid = false;
    bool fault = false;

    std::uint8_t *buffer(unsigned long int size)
    {
        if (_buffer == nullptr)
        {
            _buffer = new std::uint8_t[size];
            valid = false;
        }
        return _buffer;
    }

    ~Meta()
    {
        delete[] _buffer;
    }
};

struct Key
{
    friend class Blob;
    friend class Sample;
    friend class Remote;

    template <class ClientSample>
    friend class SamplePursuer;

    template <class ClientSample>
    friend class Request;

protected:
    union
    { // 256 bit key
        std::uint8_t bytes[32];
        std::uint_fast64_t random[4];
    };

    static std::unordered_map<Key, std::weak_ptr<Meta>, Key> metadata;

    static void vacuum() noexcept;

    std::shared_ptr<Meta> meta() const noexcept;

    const Key &operator=(const Key &that) noexcept
    {
        std::memcpy(bytes, that.bytes, 32);
        return *this;
    };

public:
    std::size_t operator()(Key const &key) const noexcept
    {                                  // relies on the fact that the
        return *(std::size_t *)(&key); // key is already a good hash
    }

    bool operator==(
        const Key &other) const noexcept
    {
        return std::memcmp(bytes, other.bytes, 32) == 0;
    }

    Key &operator=(Key &&that) noexcept
    {
        if (!(*this == that))
            std::memcpy(bytes, that.bytes, 32);
        return *this;
    }

    Key &operator=(Key &that) noexcept
    {
        if (!(*this == that))
            std::memcpy(bytes, that.bytes, 32);
        return *this;
    }

    Key(const Key &key)
    {
        std::move(key.bytes, key.bytes + 32, bytes);
    };

    Key(std::mt19937_64 &generate)
    {
        for (auto &random_bits : random)
            random_bits = generate();
    };

    Key()
    {
        std::move(no_key, no_key + 32, bytes);
    };

    explicit operator json() const noexcept;

    Key &operator=(const json &encoded);

    explicit Key(const json &encoded)
    {
        *this = encoded;
    };

    Key(std::istream &stream);
};

const Key NoKey;

class Blob
{
    friend class Remote;

    // private:
    //     static std::mutex mutex;
    //     static std::condition_variable wakeup;

protected:
    Key key;
    // size_t size();
    // std::shared_ptr<Meta> meta;

    // std::istream istream();

    // std::ostream ostream();

    // class Handler
    // {
    // public:
    //     void startFault(
    //         Key key){};
    // };

    // static Handler *handler;

public:
    // bool valid();

    // bool fault();

    // void wait();

protected:
    // static void content(
    //     Meta *meta);

    Blob(
        const Key &key);

    Blob();
};

#endif /* _EPP_METADATA_H */
