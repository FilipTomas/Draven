#include "read_preprocessing.h"


#include <memory>
#include <vector>


namespace raven{

std::uint32_t kDiamerSize = 2;
std::uint32_t kHomopolymerSize = 3;

biosoup::NucleicAcid CompressHomopolymers(biosoup::NucleicAcid &read){
    std::string sequence = read.InflateData();
    const char* seq_ptr = sequence.c_str();
    seq_ptr += 1;
    std::uint32_t repeat_count = 1;
    while(*seq_ptr != '\0'){
        if(*seq_ptr == *(seq_ptr - 1)){
            repeat_count++;
        }else{
            if(repeat_count >= kHomopolymerSize){
                sequence.erase(*(seq_ptr - repeat_count + 1), repeat_count - 1);
            }
            repeat_count = 1;
        }
    }
    return biosoup::NucleicAcid(read.name, sequence);
}

void Extended_reads::compress_homopolymers(std::shared_ptr<thread_pool::ThreadPool> thread_pool){
    // std::vector<std::future<void>> thread_futures;
    // for(std::uint32_t i; i < original_reads_.size(); i++){
    //     thread_futures.push_back(thread_pool->Submit([&](std::uint32_t j){
    //         //CompressHomopolymers(original_reads_[j]);
    //         repeat_compressed_reads_[j] = CompressHomopolymers(original_reads_[j]);
    //     }, i));
    // }
    // for(auto& it : thread_futures){
    //     it.wait();
    // }
    for (std::uint32_t i = 0; i < original_reads_.size(); i++){
        repeat_compressed_reads_.push_back(CompressHomopolymers(get_original_read(i)));
    }

}

// void Extended_reads::compress_diamers(std::shared_ptr<thread_pool::ThreadPool> thread_pool){
//     std::vector<std::future<void>> thread_futures;
//     for(std::uint32_t i; i < repeat_compressed_reads_.size(); i++){
//         thread_futures.push_back(thread_pool->Submit([&](std::uint32_t j){
//             //CompressHomopolymers(original_reads_[j]);
//             repeat_compressed_reads_[j] = CompressDiamers(repeat_compressed_reads_[j]);
//         }, i));
//     }
//     for(auto& it : thread_futures){
//         it.wait();
//     }

// }



biosoup::NucleicAcid CompressDiamers(biosoup::NucleicAcid &read){
    std::string sequence = read.InflateData();
    std::uint32_t repeat_count = 1;
    for (std::uint32_t i = 3; i < read.inflated_len; i+=2){
        if(sequence[i] == sequence[i - 2] && sequence[i - 1] == sequence[i - 3]){
            repeat_count++;
        }else{
            if(repeat_count >= kDiamerSize){
                sequence.erase(i - repeat_count*2 + 1, repeat_count*2 - 1);
                i -= repeat_count*2 - 1;
            }
            repeat_count = 1;
        }
    }
    return biosoup::NucleicAcid(read.name, sequence);
}

}