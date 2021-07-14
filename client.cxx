#include "client.h"

namespace EPP
{

    std::size_t KeyHash::operator()(Key const &key) const noexcept
    {                                  // relies on the fact that the
        return *(std::size_t *)(&key); // key is already a good hash
    }

    std::unordered_map<Key, std::weak_ptr<Meta>, KeyHash> Key::metadata;

    void Key::vacuum()
    {
        for (auto it = metadata.begin(); it != metadata.end(); it++)
            if (it->second.expired())
                metadata.erase(it->first);
    }

    std::shared_ptr<Meta> Key::meta()
    {
        static unsigned int stale = 0;
        if (stale++ > 100)
        {
            vacuum();
            stale = 0;
        }

        std::shared_ptr<Meta> strong;
        std::weak_ptr<Meta> weak;

        auto it = metadata.find(*this);
        if (it != metadata.end())
        {
            weak = it->second;
            if (weak.expired())
            {
                weak.reset();
                metadata.erase(*this);
                strong = std::shared_ptr<Meta>(new Meta());
                weak = strong;
                metadata.insert(std::pair<Key, std::weak_ptr<Meta>>(*this, weak));
            }
            else
                strong = weak.lock();
        }
        else
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
        {
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
        Parameters parameters,
        Pursuer *pursuer) noexcept
        : begin(std::chrono::steady_clock::now()), pursuer(pursuer)
    {
        _result = std::shared_ptr<Result>(new Result(parameters));
        Result *rp = _result.get();
        for (auto &random_bits : rp->key.random)
            random_bits = generate();
        pursuer->start(this);
    }

    void Server::send(const json &encoded)
    {
        std::unique_lock<std::mutex> lock(mutex);
        outgoing.push(encoded);
        wakeup.notify_one();
    }

    void Server::transmit()
    {
        std::unique_lock<std::mutex> lock(mutex);
        while (true)
            if (outgoing.empty())
                wakeup.wait(lock);
            else
            {
                json encoded = outgoing.front();
                outgoing.pop();
                Service service = request; // from json
                switch (service)
                {
                case request:
                case result:
                case fault:
                    // serialize encoded and send on the control channel
                    break;

                case content:
                    // send the blob on the data channel using stream in meta
                    break;
                }
            }
    };

    void Server::receive()
    {
        while (true)
        {
            json encoded;              // deserialize from control channel
            Service service = request; // from json
            switch (service)
            {
            case request:
                // pursuer->start(encoded);
                break;

            case result:
                // pursuer->finish(encoded);
                break;

            case fault:
                // send content service back
                send(encoded);
                break;
            }
        }
    };

    Server::Server()
    {
        transmitter = std::thread(
            [this]()
            { transmit(); });
        receiver = std::thread(
            [this]()
            { receive(); });
    };

    Server::~Server()
    {
        transmitter.join();
        receiver.join();
    }
}