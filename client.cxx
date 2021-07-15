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

    std::mutex Blob::mutex;

    std::condition_variable Blob::wakeup;

    void Blob::wait()
    {
        std::unique_lock<std::mutex> lock(mutex);
        while (!_key.meta()->valid)
            wakeup.wait(lock);
    };

    Key Sample::key()
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

    std::mt19937_64 Request::generate(random());

    void Request::finish() noexcept
    {
        pursuer->finish(this);
    }

    std::shared_ptr<Result> Request::result()
    {
        wait();
        end = std::chrono::steady_clock::now();
        _result->milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        return _result;
    }

    Request::Request(
        Pursuer *pursuer,
        Parameters parameters) noexcept
        : begin(std::chrono::steady_clock::now()), pursuer(pursuer)
    {
        _result = std::shared_ptr<Result>(new Result(parameters));
        Result *rp = _result.get();
        for (auto &random_bits : rp->key.random)
            random_bits = generate();
        pursuer->start(this);
    }

    void Pursuer::finish(
        const json &encoded)
    {
        Key request_key; // from JSON
        Request *request = requests.find(request_key)->second;
        Result *result = request->result().get();
        *result = encoded;
        request->finish();
    }
}