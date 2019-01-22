#include "util.hpp"
#include <algorithm>

void error(const std::string e) {
	std::cerr << "Error: " << e << std::endl << std::endl;
}

inline bool isalpha(const char & c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}


void strip(std::string & s) {
		for(auto it = s.begin(); it!=s.end(); ++it) {
			if(!isalpha(*it)) {
				s.erase(it);
				it--;
			}
		}
}

void write_database(const std::string & filename, const std::vector<Kmer> & initial_kmers, pCountMap * kmer2countmap, boophf_t * bphf, const ExpId2Name & exp_id2name, const ExpId2Desc & exp_id2desc, const ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Writing k-mer database to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }

	// write header
	struct HeaderDbFile hdr;
	os.write(reinterpret_cast<const char *>(&hdr.magic),sizeof(hdr.magic));
	os.write(reinterpret_cast<const char *>(&hdr.dbVer),sizeof(hdr.dbVer));

	struct HeaderDbKmers hdr_k;
	hdr_k.numKmer = initial_kmers.size();
	os.write(reinterpret_cast<const char *>(&hdr_k.numKmer),sizeof(hdr_k.numKmer));

	for(Kmer it : initial_kmers) {
		KmerIndex index = bphf->lookup(it);
		assert(index < bphf->nbKeys());
		Kmer kmer = it;
		// write k-mer
		os.write(reinterpret_cast<const char *>(&kmer),sizeof(kmer));
		// write number of experiments having this k-mer
		ExperimentCount num_exp = kmer2countmap[index]==nullptr ? 0 : static_cast<ExperimentCount>(kmer2countmap[index]->size());
		os.write(reinterpret_cast<const char *>(&num_exp),sizeof(num_exp));
		if(num_exp > 0) {
			for(auto const & it_exp : *kmer2countmap[index]) {
				ExperimentId exp_id = it_exp.first;
				assert(exp_id2name.find(exp_id) != exp_id2name.end());
				KmerCount count = kmer2countmap[index]->at(it_exp.first);
				os.write(reinterpret_cast<const char *>(&exp_id),sizeof(exp_id));
				os.write(reinterpret_cast<const char *>(&count),sizeof(count));
			}
		}
	}

	struct HeaderDbMetadata hdr_m;
	hdr_m.numExp = exp_id2name.size();
	os.write(reinterpret_cast<const char *>(&hdr_m.label),sizeof(hdr_m.label));
	os.write(reinterpret_cast<const char *>(&hdr_m.numExp),sizeof(hdr_m.numExp));

	for(auto const & it : exp_id2name) {
		const ExperimentId exp_id = it.first;
		const std::string exp_name = it.second;
		auto const it_desc = exp_id2desc.find(exp_id);
		const std::string exp_desc = ( it_desc != exp_id2desc.end() ) ? it_desc->second : "";
		assert(exp_id2readcount.count(exp_id) > 0);
		const ReadCount readcount = exp_id2readcount.at(exp_id);

		os.write(reinterpret_cast<const char *>(&exp_id),sizeof(ExperimentId));
		os.write(reinterpret_cast<const char *>(&readcount),sizeof(ReadCount));
		os.write(exp_name.c_str(),exp_name.length() + 1);
		os.write(exp_desc.c_str(),exp_desc.length() + 1);
	}

	os.close();
	if(!os) { // writing failed at some point
		error("Writing to file " + filename + "failed."); exit(EXIT_FAILURE); 
	}

}

void read_database(const std::string & filename,
										std::vector<Kmer> & initial_kmers,
										pCountMap * kmer2countmap,
										boophf_t * bphf,
										bool append,
										ExpId2Name & exp_id2name,
										ExpId2Desc & exp_id2desc,
										ExpName2Id & exp_name2id,
										ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Reading database file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::in | std::ios::binary);
	if(!ifs) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }

	//read header
	struct HeaderDbFile h_in;
	struct HeaderDbFile h_ref;
	ifs.read(reinterpret_cast<char*>(&h_in.magic), sizeof(h_in.magic));
	if(!ifs.good()) throw std::runtime_error("could not read magic bytes, file truncated");
	if(memcmp(h_in.magic,h_ref.magic,3)!=0) throw std::runtime_error("wrong file type detected");
	if(h_in.magic[3] != h_ref.magic[3]) throw std::runtime_error("file corruption detected");

	ifs.read(reinterpret_cast<char*>(&h_in.dbVer), sizeof(h_in.dbVer));
	if(!ifs.good())  throw std::runtime_error("could not read version, file truncated");

	// read k-mer section
	struct HeaderDbKmers k;
	ifs.read(reinterpret_cast<char*>(&k.numKmer), sizeof(k.numKmer));
	if(!ifs.good()) throw std::runtime_error("could not read number of kmers, file truncated");
	if(k.numKmer != bphf->nbKeys()) throw std::runtime_error("Mismatching number of k-mers in hash index and k-mer database");

	for(uint64_t n = 1; n <= k.numKmer; n++) {
		Kmer kmer;
		ifs.read(reinterpret_cast<char*>(&kmer), sizeof(Kmer));
		if(!ifs.good()) throw std::runtime_error("could not read k-mer #"+std::to_string(n)+", file truncated");
		KmerIndex index = bphf->lookup(kmer);
		assert(index < bphf->nbKeys());
		initial_kmers.emplace_back(kmer);
		ExperimentCount num_exp = 0;
		ifs.read(reinterpret_cast<char*>(&num_exp), sizeof(ExperimentCount));
		if(!ifs.good()) throw std::runtime_error("could not read experiment count for k-mer "+std::to_string(kmer)+", file truncated");
		if(num_exp>0) {
			if(append) kmer2countmap[index] = new CountMap();
			for(int i=0; i < static_cast<int>(num_exp); i++) {
				ExperimentId exp_id = 0;
				ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
				if(!ifs.good()) throw std::runtime_error("could not read experiment id for k-mer "+std::to_string(kmer)+", file truncated");
				KmerCount count = 0;
				ifs.read(reinterpret_cast<char*>(&count), sizeof(KmerCount));
				if(!ifs.good()) throw std::runtime_error("could not read k-mer count for experiment id "+std::to_string(exp_id)+", file truncated");
				if(append) kmer2countmap[index]->emplace(exp_id,count);
			}
		}
	}

	// read metadata section
	struct HeaderDbMetadata m_in;
	struct HeaderDbMetadata m_ref;

	ifs.read(reinterpret_cast<char*>(&m_in.label), sizeof(m_in.label));
	if(!ifs.good()) throw std::runtime_error("could not read metadata header, file truncated");
	if(memcmp(m_in.label,m_ref.label,8)!=0) throw std::runtime_error("invalid metadata header, file corruption detected");

	ifs.read(reinterpret_cast<char*>(&m_in.numExp), sizeof(m_in.numExp));
	if(!ifs.good()) throw std::runtime_error("could not read number of experiments in metadata section, file truncated");

	// continue reading metadata, add field for sample description
	for(uint64_t n = 1; n <= m_in.numExp; n++) {
		ExperimentId exp_id = 0;
		ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(exp_id));
		if(!ifs.good()) throw std::runtime_error("could not read experiment id for experiment #"+std::to_string(n)+", file truncated");
		assert(exp_id > 0);
		ReadCount read_count = 0;
		ifs.read(reinterpret_cast<char*>(&read_count), sizeof(read_count));
		if(!ifs.good()) throw std::runtime_error("could not read count for experiment "+std::to_string(exp_id)+", file truncated");
		std::string exp_name;
		getline(ifs, exp_name,'\0');
		if(!ifs.good()) throw std::runtime_error("could not read experiment name for experiment"+std::to_string(exp_id)+", file truncated");
		assert(exp_name.length() > 0);
		std::string exp_desc;
		getline(ifs, exp_desc,'\0');
		if(!ifs.good()) throw std::runtime_error("could not read experiment desc for experiment"+std::to_string(exp_id)+", file truncated");
		exp_id2name.emplace(exp_id,exp_name);
		exp_id2desc.emplace(exp_id,exp_desc);
		exp_id2readcount.emplace(exp_id,read_count);
		exp_name2id.emplace(exp_name,exp_id);
	}
	// there should be nothing else left after this point
	if(ifs.peek() != EOF)  throw std::runtime_error("file has extra bytes, file corruption detected");

}

void print_usage_header() {
	fprintf(stderr, "KIQ %s\n",KIQ_VERSION_STRING);
	fprintf(stderr, "Copyright 2018,2019 Peter Menzel\n");
	fprintf(stderr, "\n");
}

std::string getCurrentTime() {
	time_t t = time(0);
	char buffer[9] = {0};
	strftime(buffer, 9, "%H:%M:%S", localtime(&t));
	return std::string(buffer);
}

/**
* Converts a string of "ATCG" to a uint64_t
* where each character is represented by using only two bits
* from https://github.com/splatlab/squeakr/blob/master/kmer.h
*/
Kmer str_to_int(const std::string & str) {
	uint64_t strint = 0;
	for(auto it : str) {
		strint <<= 2;
		uint8_t curr = DNA_MAP::A;
		switch(it) {
			//case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 't': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'c': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
			case 'g': { curr = DNA_MAP::G; break; }
			// TODO right now all non-CGT chars are used as A
		}
		strint = strint | curr;
	}
	return (Kmer)strint;
}

/**
* Converts a uint64_t to a string of "ACTG"
* where each character is represented by using only two bits
* from https://github.com/splatlab/squeakr/blob/master/kmer.h
*/
std::string int_to_str(const Kmer kmer) {
	uint8_t base;
	std::string str;
	for(int i=KMER_K; i>0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = 'A';
		switch(base) {
			//	case DNA_MAP::A: { return 'A'; }
				case DNA_MAP::T: { chr = 'T'; break; }
				case DNA_MAP::C: { chr = 'C'; break; }
				case DNA_MAP::G: { chr = 'G'; break; }
		}
		str.push_back(chr);
	}
	return str;
}


void read_kmer_database(const std::string & filename, std::vector<Kmer> & initial_kmers, pCountMap * kmer2countmap, boophf_t * bphf, bool append) {
	std::cerr << getCurrentTime() << " Reading database from file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::binary);
	if(!ifs.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	uint64_t count = 0;
	while(ifs.peek()!=EOF) {
		Kmer kmer;
		ifs.read(reinterpret_cast<char*>(&kmer), sizeof(Kmer));
		KmerIndex index = bphf->lookup(kmer);
		//std::cerr << "Read k-mer " << kmer << "="<< int_to_str(kmer)  <<" with index " << index << "\n";
		assert(index < bphf->nbKeys());
		initial_kmers.emplace_back(kmer);
		ExperimentCount num_exp = 0;
		//std::cerr << "kmer index=" << index <<"\n";
		ifs.read(reinterpret_cast<char*>(&num_exp), sizeof(ExperimentCount));
		//std::cerr << "num exp=" << num_exp <<"\n";
		if(num_exp>0) {
			if(append) kmer2countmap[index] = new CountMap();
			for(int i=0; i < static_cast<int>(num_exp); i++) {
				ExperimentId exp_id = 0;
				ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
				KmerCount count = 0;
				ifs.read(reinterpret_cast<char*>(&count), sizeof(KmerCount));
				//std::cerr << "kmer_index i=" << index << " expid=" << exp_id << " count=" << count << "\n";
				if(append) kmer2countmap[index]->emplace(exp_id,count);
			}
		}
		++count;
	}
	ifs.close();
	assert(count == bphf->nbKeys());
}

void write_initial_database(const std::string & filename, const std::vector<Kmer> & initial_kmers) {
	std::cerr << getCurrentTime() << " Writing k-mer database to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }

	// write header
	struct HeaderDbFile hdr;
	os.write(reinterpret_cast<const char *>(&hdr.magic),sizeof(hdr.magic));
	os.write(reinterpret_cast<const char *>(&hdr.dbVer),sizeof(hdr.dbVer));

	struct HeaderDbKmers hdr_k;
	hdr_k.numKmer = initial_kmers.size();
	os.write(reinterpret_cast<const char *>(&hdr_k.numKmer),sizeof(hdr_k.numKmer));

	KmerCount count = 0;
	for(Kmer it : initial_kmers) {
		// write k-mer
		Kmer kmer = it;
		os.write(reinterpret_cast<const char *>(&kmer),sizeof(kmer));
		// write 0
		os.write(reinterpret_cast<const char *>(&count),sizeof(count));
	}

	struct HeaderDbMetadata hdr_m;
	hdr_m.numExp = 0;
	os.write(reinterpret_cast<const char *>(&hdr_m.label),sizeof(hdr_m.label));
	os.write(reinterpret_cast<const char *>(&hdr_m.numExp),sizeof(hdr_m.numExp));
	
	os.close();
	if(!os) { // writing failed at some point
		error("Writing to file " + filename + "failed."); exit(EXIT_FAILURE); 
	}
}

void write_kmer_database(const std::string & filename,  boophf_t * bphf,  pCountMap * kmer2countmap, const std::vector<Kmer> & initial_kmers) {

	std::cerr << getCurrentTime() << " Writing k-mer database to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	for(Kmer it : initial_kmers) {
		KmerIndex index = bphf->lookup(it);
		assert(index < bphf->nbKeys());
		Kmer kmer = it;
		// write k-mer
		os.write(reinterpret_cast<const char *>(&kmer),sizeof(Kmer));
		// write number of experiments having this k-mer
		ExperimentCount num_exp = kmer2countmap[index]==nullptr ? 0 : static_cast<ExperimentCount>(kmer2countmap[index]->size());
		os.write(reinterpret_cast<const char *>(&num_exp),sizeof(ExperimentCount));

		if(num_exp > 0) {
			for(auto const & it : *kmer2countmap[index]) {
				ExperimentId exp_id = it.first;
				KmerCount count = kmer2countmap[index]->at(it.first);
				//std::cerr << "kmer_index i=" << i << " expid=" << exp_id << " count=" << count << "\n";
				os.write(reinterpret_cast<const char *>(&exp_id),sizeof(ExperimentId));
				os.write(reinterpret_cast<const char *>(&count),sizeof(KmerCount));
			}
		}
	}
	os.close();
}


void read_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpName2Id & exp_name2id, ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Reading metadata from file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::binary);
	if(!ifs.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	while(ifs.peek()!=EOF) {
		ExperimentId exp_id = 0;
		ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
		assert(!ifs.eof());
		assert(exp_id > 0);
		ReadCount readcount = 0;
		ifs.read(reinterpret_cast<char*>(&readcount), sizeof(ReadCount));
		assert(!ifs.eof());
		std::string exp_name;
		getline(ifs, exp_name,'\0');
		assert(exp_name.length() > 0);
		exp_id2name.emplace(exp_id,exp_name);
		exp_id2readcount.emplace(exp_id,readcount);
		exp_name2id.emplace(exp_name,exp_id);
	}
	ifs.close();

}

void load_index(const std::string & filename_index,  boophf_t * bphf) {
	std::cerr << getCurrentTime() << " Reading index from file " << filename_index << "\n";
	std::ifstream ifs(filename_index, std::ios::binary);
	if(!ifs.is_open()) { std::cerr << "Cannot open file " << filename_index << std::endl; exit(EXIT_FAILURE); }
	bphf->load(ifs);
}

void write_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Writing metadata to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	for(auto const & it : exp_id2name) {
		ExperimentId exp_id = it.first;
		std::string exp_name = it.second;
		assert(exp_id2readcount.count(exp_id) > 0);
		ReadCount readcount = exp_id2readcount.at(exp_id);

		os.write(reinterpret_cast<const char *>(&exp_id),sizeof(ExperimentId));
		os.write(reinterpret_cast<const char *>(&readcount),sizeof(ReadCount));
		os.write(exp_name.c_str(),exp_name.length() + 1);
	}
	os.close();

}

