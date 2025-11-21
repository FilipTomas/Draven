// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"
#include "graph_constructor/graph_constructor.h"
#include "graph_assembler/graph_assembler.h"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "edlib.h"  // NOLINT
#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace raven::test {

  class RavenTest : public ::testing::Test {
  public:
    void SetUp() {

    }

  };

  TEST_F(RavenTest, RemoveTips) {
    raven::Graph graph;

    raven::Graph_Assembler assembler(graph);

    std::uint32_t res = assembler.RemoveTips();
    EXPECT_EQ(0, res);
  }

}  // namespace raven::test
