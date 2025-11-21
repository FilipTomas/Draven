
#ifndef RAVEN_PARSER_H
#define RAVEN_PARSER_H

#include <memory>
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"

namespace raven {

  class Parser {
  public:
    Parser() = default;

    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(const std::string &path);
  };

} // raven

#endif //RAVEN_PARSER_H
