#include "biosoup/nucleic_acid.hpp"
#include <iostream>
#include "thread_pool/thread_pool.hpp"
#include <memory>
#include <vector>


namespace raven {
    class Extended_reads{
    public:
        Extended_reads() = default;
        
        Extended_reads(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &reads, std::shared_ptr<thread_pool::ThreadPool> thread_pool) : original_reads_(std::move(reads)){
            compress_homopolymers(thread_pool);
//            compress_diamers(thread_pool);
        };

        std::vector<biosoup::NucleicAcid> get_compressed_reads(){
            return repeat_compressed_reads_;
        }

        // biosoup::NucleicAcid get_compressed_read(std::uint32_t i){
        //     return *repeat_compressed_reads_[i];
        // }

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> get_original_reads(){
            return original_reads_;
        }

        biosoup::NucleicAcid get_original_read(std::uint32_t i){
            return *original_reads_[i];
        }


    private:
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> original_reads_;
    std::vector<biosoup::NucleicAcid> repeat_compressed_reads_{original_reads_.size()};
    std::vector<std::uint32_t> homopolymer_positions_;
    std::vector<std::uint32_t> diamer_positions_;
    void compress_diamers(std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    void compress_homopolymers(std::shared_ptr<thread_pool::ThreadPool> thread_pool);
    };

}