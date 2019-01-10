#pragma once

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include <time.h>
#include <map>
#include <fstream>

#include "BooPHF/BooPHF.h"
#include "version.hpp"

using hasher_t = boomphf::SingleHashFunctor<u_int64_t>;
using boophf_t = boomphf::mphf<u_int64_t, hasher_t>;

using ExperimentId = uint32_t;
using ExperimentCount = uint32_t;
using KmerCount = uint32_t;
using KmerIndex = uint64_t;
using Kmer = uint64_t;
using ReadCount = uint64_t;

using CountMap = std::map<ExperimentId, KmerCount>;
using pCountMap = CountMap*;
using ExpId2Name = std::map<ExperimentId, std::string>;
using ExpName2Id = std::map<std::string, ExperimentId>;
using ExpId2ReadCount = std::map<ExperimentId, ReadCount>;

enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

#define KMER_K 32

void error(const std::string e);

void strip(std::string &s);

bool isalpha(const char & c);

void print_usage_header();

std::string getCurrentTime();

Kmer str_to_int(const std::string & str);
std::string int_to_str(Kmer kmer);

void load_index(const std::string & filename_index,  boophf_t * bphf);

void read_kmer_database(const std::string & filename,std::vector<Kmer> & initial_kmers, pCountMap * kmer2countmap, boophf_t * bphf, bool append);
void write_kmer_database(const std::string & filename, boophf_t * bphf, pCountMap * kmer2countmap, const std::vector<Kmer> & initial_kmers);
void write_initial_kmer_database(const std::string & filename, const std::vector<Kmer> & initial_kmers);

void read_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpId2ReadCount & exp_id2readcount);
void read_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpName2Id & exp_name2id, ExpId2ReadCount & exp_id2readcount);
void write_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpId2ReadCount & exp_id2readcount);

