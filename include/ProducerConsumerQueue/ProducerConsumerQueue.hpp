/*

This implementation of a Producer-Consumer-Queue is adapted from
ProducerConsumerQueue (v0.1) by Lasse Maretty and Jonas Andreas Sibbesen

Major changes are:
- put all code in a single header file
- push() only notifies consumers once the queue is sufficiently full, which
	avoids slowdown due to constant switching between consumers and producers
	when production is slower than consumption
- more code cleanup

Copyright (c) 2018 Peter Menzel
Copyright (c) 2014 Lasse Maretty and Jonas Andreas Sibbesen


The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/


#ifndef __ProducerConsumerQueue__
#define __ProducerConsumerQueue__

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <atomic>
#include <math.h>

template<typename Data>
class ProducerConsumerQueue {

	private:
		size_t max_buffer_size;
		size_t notify_buffer_size;
		static constexpr float notify_fill_factor = 4.0 / 5.0;
		bool pushed_last = false;
		std::deque<Data> queue;
		std::mutex queue_mutex;
		std::condition_variable producer_cv;
		std::condition_variable consumer_cv;

	public:
		ProducerConsumerQueue(const ProducerConsumerQueue&) = delete;
		ProducerConsumerQueue& operator=(const ProducerConsumerQueue&) = delete;
		ProducerConsumerQueue(size_t max_buffer_size_) : max_buffer_size(max_buffer_size_)  {
			notify_buffer_size = floor(max_buffer_size * notify_fill_factor);
		}

		void push(const Data & new_element) {
			std::unique_lock<std::mutex> queue_lock(queue_mutex);

				producer_cv.wait(queue_lock,[this]{return queue.size() < max_buffer_size;});

				queue.push_back(new_element);

			queue_lock.unlock();

			if(queue.size() > notify_buffer_size)
				consumer_cv.notify_one();
		}

		bool pop(Data & returned_element) {
			std::unique_lock<std::mutex> queue_lock(queue_mutex);

				if(pushed_last && queue.empty()) {
					return false;
				}

				consumer_cv.wait(queue_lock,[this]{return !queue.empty();});

				returned_element = queue.front();
				queue.pop_front();

			queue_lock.unlock();

			producer_cv.notify_one();

			return true;
		}

		void pushedLast() {
			std::unique_lock<std::mutex> queue_lock(queue_mutex);
				pushed_last = true;
			queue_lock.unlock();
			consumer_cv.notify_all();
		}

		size_t size() {
			std::lock_guard<std::mutex> g(queue_mutex);
			return queue.size();
		}

};

#endif
