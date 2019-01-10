#include "CountThread.hpp"
#include "util.hpp"

CountThread::CountThread(queue_t * queue_, std::atomic<KmerCount> * counts_, boophf_t * bphf_, std::unordered_set<Kmer> * initial_kmers_) noexcept {
	queue = queue_;
	counts = counts_;
	bphf = bphf_;
	n_elem = (KmerIndex)bphf->nbKeys();
	initial_kmers = initial_kmers_;
	assert(initial_kmers->size() == n_elem);
}

void CountThread::doWork() noexcept {
	ReadItem * item = nullptr;
	while(queue->pop(item)) {
		assert(item != nullptr);
		assert(item->sequence1.length() >= KMER_K);
		uint64_t t = 0;
		const char * c = item->sequence1.c_str();
		//init first K chars
		for(int k=0; k<KMER_K; k++) {
			uint8_t curr = DNA_MAP::A;
			switch (*c) {
				//case 'A': { curr = DNA_MAP::A; break; }
				case 'T': { curr = DNA_MAP::T; break; }
				case 'C': { curr = DNA_MAP::C; break; }
				case 'G': { curr = DNA_MAP::G; break; }
			}
			t = t | curr;
			if(k < KMER_K - 1) {
				t = t << 2;
				c++;
			}
		}
		while(true) {
			if(initial_kmers->find(t) != initial_kmers->end()) {
				const KmerIndex index = bphf->lookup(t);
				assert(index < n_elem);
				counts[index]++;
			}
			c++;
			if(*c == '\0') break;
			uint8_t curr = DNA_MAP::A;
			switch (*c) {
				//case 'A': { curr = DNA_MAP::A; break; }
				case 'T': { curr = DNA_MAP::T; break; }
				case 'C': { curr = DNA_MAP::C; break; }
				case 'G': { curr = DNA_MAP::G; break; }
			}
			t = t << 2;
			t = t | curr;
		}

		delete item;

	} // end while queue

}


