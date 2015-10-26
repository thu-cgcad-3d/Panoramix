#pragma once

#include <thread>

namespace pano {

    namespace core {

        template <class FunT>
        void ParallelRun(int n, int batchSize, FunT && fun) {
            std::vector<std::thread> threads;
            threads.reserve(batchSize);
            for (int i = 0; i < n; i++) {
                threads.emplace_back(std::forward<FunT>(fun), i);
                if (threads.size() >= batchSize || i == n-1) {
                    for (auto & t : threads) {
                        t.join();
                    }
                    threads.clear();
                }
            }          
        }

    }

}