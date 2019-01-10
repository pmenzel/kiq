#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
#include <array>
#include <string>
#include <atomic>
#include <deque>
#include <stdexcept>
#include <memory>

#include <ncbi-vdb/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include "ProducerConsumerQueue/ProducerConsumerQueue.hpp"
#include "BooPHF/BooPHF.h"
#include "util.hpp"
#include "ReadItem.hpp"
#include "CountThread.hpp"


void usage_ksra() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq sra -i <file> -k <file> -m <file> -l <file> [-a] [-z <int>]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of index file\n");
	fprintf(stderr, "   -k FILENAME   Name of k-mer count database file\n");
	fprintf(stderr, "   -m FILENAME   Name of sample metadata file\n");
	fprintf(stderr, "   -l FILENAME   Name of file with sample list in TSV format \n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -a            Append mode\n");
	fprintf(stderr, "   -z INT        Fixed number of parallel threads for counting (default: 5)\n");
	fprintf(stderr, "   -v            Enable verbose output\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

int main_ksra(int argc, char** argv) {

	size_t max_num_threads = 5;
	size_t curr_num_threads = 5;
	size_t max_queue_size = 9999;
	const float ma_alpha = 0.7f;
	bool append = false;
	bool debug = false;
	bool verbose = false;

	std::string filename_index;
	std::string filename_db;
	std::string filename_meta;
	std::string filename_inputlist;

	ncbi::NGS::setAppVersionString("kiq-0.1");
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvai:k:m:l:z:")) != -1) {
		switch (c)  {
			case 'h':
				usage_ksra();
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'a':
				append = true; break;
			case 'm':
				filename_meta = optarg; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'l':
				filename_inputlist = optarg; break;
			case 'z':
				max_num_threads = atoi(optarg);
				curr_num_threads = max_num_threads;
				break;
			default:
				usage_ksra();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_ksra(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_ksra(); }
	if(filename_meta.length() == 0) { error("Please specify the name of the metadata file, using the -m option."); usage_ksra(); }
	if(filename_inputlist.length() == 0) { error("Please specify the name of the sample list file, using the -l option."); usage_ksra(); }

	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>();
	load_index(filename_index, bphf);

	KmerIndex n_elem = (KmerIndex)bphf->nbKeys();
	std::vector<Kmer> initial_kmers;
	pCountMap * kmer2countmap = new pCountMap[n_elem](); // init new array of size n_elem

	read_kmer_database(filename_db,initial_kmers,kmer2countmap,bphf,append);

	std::unordered_set<Kmer> initial_kmers_set;
	initial_kmers_set.insert(initial_kmers.cbegin(),initial_kmers.cend());

	ExpId2Name exp_id2name;
	ExpName2Id exp_name2id;
	ExpId2ReadCount exp_id2readcount;

	ExperimentId experiment_numericid = 0;

	if(append) {
		read_experiment_database(filename_meta, exp_id2name, exp_name2id, exp_id2readcount);
		experiment_numericid = (ExperimentId)exp_id2name.size();
	}

	std::atomic<KmerCount> * tmp_counts_atomic = new std::atomic<KmerCount>[n_elem];

	// FOR EACH EXPERIMENT

	std::ifstream ifs_inputlist(filename_inputlist);
	if(!ifs_inputlist.is_open()) { std::cerr << "Cannot open file " << filename_inputlist << std::endl; exit(EXIT_FAILURE); }

	std::string line;
	while(getline(ifs_inputlist, line)) {
		if(line.length() == 0) { continue; }
		size_t tab = line.find('\t');
		if(tab == std::string::npos || tab == 0) {
			continue;
		}
		std::string experiment_stringid = line.substr(0,tab);
		std::string filename_seq = line.substr(tab+1);

		if(exp_name2id.count(experiment_stringid) > 0) { // experiment name is already in DB
			if(append) {
				std::cerr << "Experiment " << experiment_stringid << " is already in database. Skipping...\n" ;
				continue; // don't do that experiment again in append mode
			}
			experiment_numericid = exp_name2id.at(experiment_stringid);
		} else {
			experiment_numericid = (ExperimentId)exp_id2name.size() + 1;
			exp_id2name.emplace(experiment_numericid,experiment_stringid);
			exp_name2id.emplace(experiment_stringid,experiment_numericid);
		}

		memset(tmp_counts_atomic,0,n_elem*sizeof(std::atomic<KmerCount>));

		auto queue = new queue_t(max_queue_size);
		std::deque<std::thread> threads;
		std::deque<std::unique_ptr<CountThread>> threadpointers;
		for(int i=0; i < static_cast<int>(curr_num_threads); i++) {
			std::unique_ptr<CountThread> p(new CountThread(queue,tmp_counts_atomic,bphf,&initial_kmers_set));
			threads.push_back(std::thread(&CountThread::doWork,p.get()));
			threadpointers.push_back(std::move(p));
		}

		/* vars for progress indicator */
		unsigned char x[4] = { '/','-','\\','|'};
		size_t count=0;
		float ma = 1; // running average of queue size
		int i = 0;

		ngs::ReadCollection sra_read_coll = ncbi::NGS::openReadCollection(filename_seq.c_str());
		ngs::ReadIterator sra_it = sra_read_coll.getReadRange(1, sra_read_coll.getReadCount(), ngs::Read::all);

		std::cerr << getCurrentTime() << " Start counting k-mers from "<< experiment_stringid << "\n";

		if(sra_it.nextRead()) // go to first read
		while(sra_it.nextFragment() || (sra_it.nextRead() && sra_it.nextFragment())) { // either go to next fragment (2nd read in pair) in the current read or advance to next read and its first seqgment
			ngs::StringRef bases = sra_it.getFragmentBases();
			const char * pSeq = bases.data();
			uint64_t l = bases.size();
			std::string sequence(pSeq,l);
			//std::cout << sequence << "\n";
			strip(sequence); // remove non-alphabet chars
			if(sequence.length()>=K) {
				queue->push(new ReadItem(sequence));
			}
			if(count++ % 100000 == 0) {
				size_t curr_queue_size = queue->size();
				ma = ma * (1-ma_alpha) + static_cast<float>(curr_queue_size) * ma_alpha;
				fprintf(stderr,"\r%c %lu/%lu %4i ", x[i++], curr_num_threads, max_num_threads, static_cast<int>(ma));
				if(i>3) i=0;
			}
		}
		fprintf(stderr,"\r                              \r");
		fflush(stderr);

		// finish reading file

		queue->pushedLast();
		while(!threads.empty()) {
			threads.front().join();
			threads.pop_front();
			threadpointers.pop_front();
		}
		delete queue;

		std::cerr << getCurrentTime() << " Processed " << count << " sequences\n";

		if(exp_id2readcount.count(experiment_numericid) > 0) {
			exp_id2readcount.at(experiment_numericid) =  exp_id2readcount.at(experiment_numericid) + count;
		}
		else {
			exp_id2readcount.emplace(experiment_numericid,count);
		}

		//then finally go through all k-mers and add to big hash map
		for(KmerIndex i=0; i < n_elem; i++) {
			if(tmp_counts_atomic[i]==0) {
				continue;
			}
			if(kmer2countmap[i]==nullptr) { // experiment id has not been seen for this k-mer before
				kmer2countmap[i] = new CountMap();
				kmer2countmap[i]->emplace(experiment_numericid, tmp_counts_atomic[i]);
			}
			else {
				auto pos = kmer2countmap[i]->find(experiment_numericid);
				if(pos != kmer2countmap[i]->end()) {
					pos->second += tmp_counts_atomic[i];
				}
				else {
					kmer2countmap[i]->emplace(experiment_numericid, tmp_counts_atomic[i]);
				}
			}
		}

		// SAVE DATABASE to file
		write_kmer_database(filename_db, bphf, kmer2countmap, initial_kmers);

		// save experiment metadata to file
		write_experiment_database(filename_meta,exp_id2name,exp_id2readcount);

	} // end while list of all experiments to read from files


	delete[] tmp_counts_atomic;
	for(KmerIndex i = 0; i < n_elem;i++) {
		if(kmer2countmap[i] != nullptr) {
			delete kmer2countmap[i];
		}
	}
	delete[] kmer2countmap;
	delete bphf;

	return 0;

}


