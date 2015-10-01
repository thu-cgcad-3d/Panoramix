#pragma once

#include "../core/basic_types.hpp"

namespace pano {
    namespace misc {


        std::string Tagify(const std::string & path);

        template <class StringT, class ... Ts>
        inline bool SaveCache(const std::string & path, StringT && what, Ts && ... ts) {
            return pano::core::SaveToDisk("./cache/" + Tagify(path) + "_" + what + ".cereal", ts...);
        }

        template <class StringT, class ... Ts>
        inline bool LoadCache(const std::string & path, StringT && what, Ts & ... ts) {
            return pano::core::LoadFromDisk("./cache/" + Tagify(path) + "_" + what + ".cereal", ts...);
        }


    }
}