// Copyright (c) 2020 Robert Vaser

#include <iostream>
#include <fstream>
#include <filesystem>

#include "biosoup/timer.hpp"

#include "option_manager.h"
#include "graph.hpp"
#include "graph_constructor/graph_constructor.h"
#include "graph_assembler/graph_assembler.h"
#include "parser.h"

//#include "catboost_model.hpp"
#include "cmath"

#ifdef KMC_BIN_DIR
  #include "kmc_file.h"
  #include "kmer_api.h"
#endif
#include <sstream>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};


int main(int argc, char **argv) {
  Program_Parameters param{};

  if (!ProcessParameters(argc, argv, param))
    return 0;

  raven::min_unitig_size = param.min_unitig_size;

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(param.num_threads);

  raven::Graph graph{param.checkpoints};

  if (param.resume) {
    try {
      graph.Load(param.cereal_filename);
    } catch (std::exception &exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded previous run "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();
  }

  raven::Parser parser;
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;

  if ((graph.state_manager_.construct_any_state() || param.num_polishing_rounds > 2) && !param.skip_contruction) {
    try {
      auto sparser = parser.CreateParser(param.sequence_path);
      if (sparser == nullptr) {
        return 1;
      }
      sequences = sparser->Parse(-1);
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (sequences.empty()) {
      std::cerr << "[raven::] error: empty sequences set" << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded " << sequences.size() << " sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

    timer.Start();
  } else {
      std::cerr << "[raven::] skipped sequence loading" << std::endl;
  }

  raven::Graph_Constructor graph_constructor{graph, thread_pool};
  if (!param.skip_contruction){
    std::cout << "Constructing graph with params: kmer_size:" << param.kmer_len  << " winodw_size:" << param.window_len << " " << std::endl;
    graph_constructor.Construct(sequences, param);
  } else if (param.load_input_gfa) {
    graph_constructor.LoadFromGfa(param.input_gfa_path);
  } else if (param.load_input_paf) {
    graph_constructor.LoadFromPaf(sequences, param.input_paf_path);
  } else {
    std::cerr << "[raven::] error: unknown option" << std::endl;
    return 1;
  }

  graph.PrintGfa(param.gfa_post_construction_filename, param.print_gfa_seq);

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> ul_sequences;
  if (!param.ul_read_path.empty()){
    try {
      auto ul_sequence_parser = parser.CreateParser(param.ul_read_path);
      ul_sequences = ul_sequence_parser->Parse(-1);
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
    }

    if (ul_sequences.empty()) {
      std::cerr << "[raven::] error: ul read path set but the file appears empty\n" 
                << "[raven::] reverting to no ul assembly"<< std::endl;
    } else{
      std::cerr << "[raven::] loaded " << ul_sequences.size() << " ul sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();
    }
  }

  raven::Graph_Assembler graph_assembler{graph, param, thread_pool};
  if (ul_sequences.empty()) {
    graph_assembler.Assemble();
  } else {
    graph_assembler.UlAssemble(ul_sequences);
  }

  graph.PrintGfa(param.gfa_post_cleaning_filename, param.print_gfa_seq);
  if(param.ploidy >= 2){
    graph_assembler.AssembleDiploids();
  } else {
    graph_assembler.AssembleHaploids();
  }

  graph.PrintGfa(param.gfa_path, param.print_gfa_seq);

  if(param.ploidy >= 2) {
    // output to file
    std::filesystem::path path1 = param.root_path;
    std::filesystem::path path2 = param.root_path;
    path1 += "-1.fasta";
    path2 += "-2.fasta";
    std::ofstream outfile1, outfile2;
    outfile1.open(path1);
    if (!outfile1.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path1 << std::endl;
      return 1;
    }
    outfile2.open(path2);
    if (!outfile2.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path2 << std::endl;
      outfile1.close();
      return 1;
    }

    for (const auto &it: graph.GetAssembledData(true)) {
      outfile1 << ">" << it->name << std::endl;
      outfile1 << it->InflateData() << std::endl;
    }

    for (const auto &it: graph.GetAssembledData(false)) {
      outfile2 << ">" << it->name << std::endl;
      outfile2 << it->InflateData() << std::endl;
    }

    outfile1.close();
    outfile2.close();
  }
  else if(param.ploidy == 1) {
    std::filesystem::path path1 = param.root_path;
    path1 += ".fasta";
    std::ofstream outfile1;

    outfile1.open(path1);
    if (!outfile1.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path1 << std::endl;
      return 1;
    }

    for (const auto &it: graph.GetAssembledData(true)) {
      outfile1 << ">" << it->name << std::endl;
      outfile1 << it->InflateData() << std::endl;
    }

    outfile1.close();
  }

  timer.Stop();
  std::cerr << "[raven::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
