#include "client.h"

namespace EPP
{
    /**
     * implementation details for client interface
     */

    std::size_t KeyHash::operator()(Key const &key) const noexcept
    {                                  // relies on the fact that the
        return *(std::size_t *)(&key); // key is already a good hash
    }

    std::unordered_map<Key, std::weak_ptr<Meta>, KeyHash> Key::metadata;

    void Key::vacuum() noexcept
    {
        // remove expired links
        for (auto it = metadata.begin(); it != metadata.end(); it++)
            if (it->second.expired())
                metadata.erase(it->first);
    }

    std::shared_ptr<Meta> Key::meta() const noexcept
    {
        // periodically clean up
        static unsigned int stale = 0;
        if (stale++ > 100)
        {
            vacuum();
            stale = 0;
        }

        std::shared_ptr<Meta> strong;
        std::weak_ptr<Meta> weak;
        // look to see if it's still in memory
        auto it = metadata.find(*this);
        if (it != metadata.end())
        {
            weak = it->second;
            if (weak.expired()) // nope, expired
            {
                weak.reset();
                metadata.erase(*this);
            }
            else
                strong = weak.lock(); // yep, take ownership
        }
        // create one if necessary
        if (!strong)
        {
            strong = std::shared_ptr<Meta>(new Meta());
            weak = strong;
            metadata.insert(std::pair<Key, std::weak_ptr<Meta>>(*this, weak));
        }

        return strong;
    }

    Key::Key(std::istream &stream)
    {
        Key block;
        do
        { // really should be SHA256
            stream.get((char *)(&block), 32);
            for (int i = 0; i < 4; i++)
                random[i] ^= block.random[i];
        } while (!stream.eof());
    };

    bool Blob::valid()
    {
        std::unique_lock<std::mutex> lock(mutex);
        return meta->valid;
    };

    bool Blob::fault()
    {
        std::unique_lock<std::mutex> lock(mutex);
        if (meta->valid || meta->fault)
            return false;

        handler->startFault(_key);
        meta->fault = true;
        return true;
    };

    void Blob::wait()
    {
        std::unique_lock<std::mutex> lock(mutex);
        if (meta->valid)
            return;
        if (!meta->fault)
        {
            handler->startFault(_key);
            meta->fault = true;
        }
        while (!meta->valid)
            wakeup.wait(lock);
    };

    void Blob::content(
        Meta *meta)
    {
        std::unique_lock<std::mutex> lock(mutex);
        meta->fault = false;
        meta->valid = true;
        wakeup.notify_all();
    };

    Blob::Blob(const Key &key) : _key(key), meta(_key.meta()){};

    Blob::Blob() = default;

    std::mutex Blob::mutex;

    std::condition_variable Blob::wakeup;

    Blob::Handler *Blob::handler;

    const Key &Sample::key()
    {
        SampleStream *stream = new SampleStream(*this);
        if (_key == NoKey)
        {
            Key real_key(stream);
            _key = real_key;
            stream->seekg(0);
        }
        // meta().stream = stream;

        return _key;
    }

    Key Subset::key()
    {
        SubsetStream *stream = new SubsetStream(*this);
        if (_key == NoKey)
        {
            Key real_key(stream);
            _key = real_key;
            stream->seekg(0);
        }
        // meta().stream = stream;

        return _key;
    }

    std::mt19937_64 _Request::generate(random());

    Request::Request(
        Pursuer *pursuer,
        Parameters parameters)
    {
        *this = pursuer->start(new _Request(pursuer, parameters));
    };

    const Key Request::key() const noexcept
    {
        return (*this)->working_result->key;
    };

    void Request::finish()
    {
        (*this)->pursuer->finish(get());
    }

    bool Request::finished() const noexcept
    {
        return (*this)->_finished;
    };

    void Request::wait() const noexcept
    {
        (*this)->pursuer->wait(get());
    };

    Result Request::result() const noexcept
    {
        if (!finished())
            wait();
        if (!(*this)->final_result)
        {
            (*this)->end = std::chrono::steady_clock::now();
            (*this)->working_result->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>((*this)->end - (*this)->begin);
            (*this)->final_result = std::shared_ptr<_Result>((*this)->working_result);
        }
        return (*this)->final_result;
    }

    Request& Request::operator++()
    {
        std::unique_lock<std::mutex> lock((*this)->pursuer->mutex);
        ++(*this)->outstanding;
        return *this;
    };

    Request& Request::operator--()
    {
        {
            std::unique_lock<std::mutex> lock((*this)->pursuer->mutex);
            --(*this)->outstanding;
        }
        if ((*this)->outstanding == 0)
            finish();
        return *this;
    };

    _Request::_Request(
        Pursuer *pursuer,
        Parameters parameters) noexcept
        : begin(std::chrono::steady_clock::now()), pursuer(pursuer),
        _finished(false), outstanding(0)
    {
        working_result = new _Result(parameters);
        for (auto &random_bits : working_result->key.random)
            random_bits = generate();
    }

    Request Pursuer::start(
        _Request *request) noexcept
    {
        Request shared(request);
        std::unique_lock<std::mutex> lock(mutex);
        bool inserted = requests.insert(std::pair<const Key, Request>(request->working_result->key, shared)).second;
        assert(inserted);
        request->_finished = false;
        return shared;
    }

    void Pursuer::finish(
        _Request *request) noexcept
    {
        std::unique_lock<std::mutex> lock(mutex);
        auto it = requests.find(request->working_result->key);
        assert(it != requests.end());
        requests.erase(it);
        request->_finished = true;
        request->completed.notify_all();
    };

    void Pursuer::wait(
        _Request *request) noexcept
    {
        std::unique_lock<std::mutex> lock(mutex);
        while (!request->_finished)
            request->completed.wait(lock);
    };

    void Pursuer::finish(
        const json &encoded)
    {
        Key request_key; // from JSON
        Request request = requests.find(request_key)->second;
        _Result *result = request.working();
        *result = encoded;
        request.finish();
    }
}