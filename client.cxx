#include "client.h"

namespace EPP
{
    /**
     * implementation details for client interface
     */


    // std::unordered_map<Key, std::weak_ptr<Meta>, KeyHash> Key::metadata;

    // void Key::vacuum() noexcept
    // {
    //     // remove expired links
    //     for (auto it = metadata.begin(); it != metadata.end(); it++)
    //         if (it->second.expired())
    //             metadata.erase(it->first);
    // }

    // std::shared_ptr<Meta> Key::meta() const noexcept
    // {
    //     // periodically clean up
    //     static unsigned int stale = 0;
    //     if (stale++ > 100)
    //     {
    //         vacuum();
    //         stale = 0;
    //     }

    //     std::shared_ptr<Meta> strong;
    //     std::weak_ptr<Meta> weak;
    //     // look to see if it's still in memory
    //     auto it = metadata.find(*this);
    //     if (it != metadata.end())
    //     {
    //         weak = it->second;
    //         if (weak.expired()) // nope, expired
    //         {
    //             weak.reset();
    //             metadata.erase(*this);
    //         }
    //         else
    //             strong = weak.lock(); // yep, take ownership
    //     }
    //     // create one if necessary
    //     if (!strong)
    //     {
    //         strong = std::shared_ptr<Meta>(new Meta());
    //         weak = strong;
    //         metadata.insert(std::pair<Key, std::weak_ptr<Meta>>(*this, weak));
    //     }

    //     return strong;
    // }

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

    // bool Blob::valid()
    // {
    //     std::unique_lock<std::mutex> lock(mutex);
    //     return meta->valid;
    // };

    // bool Blob::fault()
    // {
    //     std::unique_lock<std::mutex> lock(mutex);
    //     if (meta->valid || meta->fault)
    //         return false;

    //     handler->startFault(_key);
    //     meta->fault = true;
    //     return true;
    // };

    // void Blob::wait()
    // {
    //     std::unique_lock<std::mutex> lock(mutex);
    //     if (meta->valid)
    //         return;
    //     if (!meta->fault)
    //     {
    //         handler->startFault(_key);
    //         meta->fault = true;
    //     }
    //     while (!meta->valid)
    //         wakeup.wait(lock);
    // };

    // void Blob::content(
    //     Meta *meta)
    // {
    //     std::unique_lock<std::mutex> lock(mutex);
    //     meta->fault = false;
    //     meta->valid = true;
    //     wakeup.notify_all();
    // };

    // Blob::Blob(const Key &key) : _key(key), meta(_key.meta()){};

    Blob::Blob() = default;

    // std::mutex Blob::mutex;

    // std::condition_variable Blob::wakeup;

    // Blob::Handler *Blob::handler;

    // const Key &Sample::key() noexcept
    // {
    //     SampleStream *stream = new SampleStream(*this);
    //     if (_key == NoKey)
    //     {
    //         Key real_key(*stream);
    //         // _key = real_key;
    //         stream->seekg(0);
    //     }
    //     // meta().stream = stream;

    //     return _key;
    // }

    // Request::Request(
    //     Pursuer *pursuer,
    //     Parameters parameters){
    //     // *this = pursuer->start(new _Request(pursuer, parameters));
    // };

    // bool Request::finished() const noexcept
    // {
    //     return (*this)->_finished;
    // };

    // void Request::wait() const noexcept
    // {
    //     (*this)->pursuer->wait(get());
    // };

    // Result Request::result() const noexcept
    // {
    //     if (!finished())
    //         wait();

    //     return (*this)->result;
    // }

    const unsigned short Parameters::N = 1 << 8; // resolution of points and boundaries

    void Remote::out(const json &encoded) // does not block
    {
        std::unique_lock<std::mutex> lock(mutex);
        outgoing.push(encoded);
        wake_out.notify_one();
    }

    json Remote::in() // blocks calling thread
    {
        std::unique_lock<std::mutex> lock(mutex);
        while (incoming.empty())
            wake_in.wait(lock);
        json encoded = incoming.front();
        incoming.pop();
        return encoded;
    }

    void Remote::copy(
        std::istream *in,
        std::ostream *out,
        unsigned long int count)
    {
        char buffer[8192];
        while (count > 0)
        {
            unsigned long int chunk = count;
            if (chunk > 8192)
                chunk = 8192;
            in->read(buffer, chunk);
            long int n = in->gcount();
            if (n > 0)
                out->write(buffer, n);
            count -= n;
        }
    };

    void Remote::startFault(
        Key key)
    {
        json encoded;
        out(encoded);
    };

    void Remote::transmit()
    {
        // don't need to hold the lock while we do I/O
        json encoded;
        {
            std::unique_lock<std::mutex> lock(mutex);
            if (outgoing.empty())
            {
                wake_out.wait(lock);
                return; // make sure we're still on the air
            }
            else
            {
                encoded = outgoing.front();
                outgoing.pop();
            }
        }
        Service service = request; // from json
        switch (service)
        {
        case request:
        case result:
        case fault:
        {
            // serialize encoded and send on the control channel
            std::string serialized("<our json message>");
            remote_control->write(serialized.data(), serialized.size());
            break;
        }
        case content:
        {
            // send the blob on the data channel using stream in meta
            Key blob_key; // from json
            // Meta *meta = blob_key.meta().get();
            // std::iostream *blob = meta->stream;
            // blob->clear();
            // blob->seekp(0);

            // serialize encoded and send on the control channel
            std::string serialized("Whatever");
            remote_control->write(serialized.data(), serialized.size());
            // send the blob content on the data channel
            // copy(blob, remote_data, meta->size);
            remote_data->flush();
            break;
        }
        }
        remote_control->flush();
    };

    void Remote::receive()
    {
        json encoded;              // deserialize from control channel
        Service service = request; // from json
        switch (service)
        {
        case request:
        case result:
        { // kick upstairs
            std::unique_lock<std::mutex> lock(mutex);
            incoming.push(encoded);
            wake_in.notify_one();
        }
        case content:
        {
            Key blob_key; // from json
            // Meta *meta = blob_key.meta().get();
            // std::ostream *content = meta->stream;
            // content->clear();
            // content->seekp(0);
            // copy(remote_data, content, meta->size);
            // content->flush();
            // Blob::content(meta); // wakeup anyone waiting for the data
            break;
        }
        case fault:
        { // change json from fault to content service and
            // put it on the output queue
            out(encoded);
            break;
        }
        }
    };

    Remote::Remote()
    {
        // Blob::handler = this;
        transmitter = std::thread(
            [this]()
            {
                while (on_the_air)
                    transmit();
            });
        receiver = std::thread(
            [this]()
            {
                while (on_the_air)
                    receive();
            });
    };

    Remote::~Remote()
    {
        {
            std::unique_lock<std::mutex> lock(mutex);
            on_the_air = false;
            wake_in.notify_all();
            wake_out.notify_all();
        }
        transmitter.join();
        receiver.join();
    };
}