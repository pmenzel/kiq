#pragma once

#include <stdint.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <atomic>
#include <cstring>
#include <climits>
#include <map>
#include <utility>
#include <functional>
#include <locale>

#include "ProducerConsumerQueue/ProducerConsumerQueue.hpp"
#include "ReadItem.hpp"
#include "util.hpp"

using queue_t = ProducerConsumerQueue<ReadItem*>;

class CountThread {
	protected:
	queue_t * queue;
	std::atomic<KmerCount> * counts;
	boophf_t * bphf;
	KmerIndex n_elem;
	std::unordered_set<Kmer> * initial_kmers;

	public:
	CountThread(queue_t *, std::atomic<KmerCount> *,boophf_t *, std::unordered_set<Kmer> *) noexcept;
	void doWork() noexcept;
	CountThread(CountThread const&) = delete;
	void operator=(CountThread const&) = delete;


};

