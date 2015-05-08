#ifndef PANOLYZ_ROUTINES_HPP
#define PANOLYZ_ROUTINES_HPP

#include "../../src/core/basic_types.hpp"

namespace panolyz {

#define DECL_ALGO(T) namespace T {void Run();}

    // declare routines
    DECL_ALGO(PanoramaIndoor);
    //DECL_ALGO(YorkUrbanDB);
    //DECL_ALGO(YorkUrbanDB2);
    DECL_ALGO(NormalIndoor);
    DECL_ALGO(NYU2);
    DECL_ALGO(Prepare);

    inline std::string Tagify(const std::string & path){
        auto tag = path;
        std::replace(tag.begin(), tag.end(), '.', '_');
        std::replace(tag.begin(), tag.end(), '\\', '.');
        std::replace(tag.begin(), tag.end(), '/', '.');
        std::replace(tag.begin(), tag.end(), ':', '.');
        return tag;
    }


    template <class StringT, class ... Ts>
    inline void Save(const std::string & path, StringT && s, Ts && ... ts){
        panoramix::core::SaveToDisk("./cache/" + Tagify(path) + "_" + s + ".cereal", ts...);
    }

    template <class StringT, class ... Ts>
    inline void Load(const std::string & path, StringT && s, Ts & ... ts){
        panoramix::core::LoadFromDisk("./cache/" + Tagify(path) + "_" + s + ".cereal", ts...);
    }




}

 
#endif