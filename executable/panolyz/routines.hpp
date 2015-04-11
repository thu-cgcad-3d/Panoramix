#ifndef PANOLYZ_ROUTINES_HPP
#define PANOLYZ_ROUTINES_HPP

#include "../../src/core/basic_types.hpp"

namespace panolyz {

    enum PanolyzAlgorithm {
        PanoramaIndoor,
        YorkUrbanDB,
        PanoContext
    };

    template <PanolyzAlgorithm Algo>
    void Routine(std::integral_constant<PanolyzAlgorithm, Algo> = {});

#define ROUTINE_FOR_ALGORITHM(Algo) \
    template <> \
    void Routine(std::integral_constant<PanolyzAlgorithm, Algo>)


    // declare routines

    ROUTINE_FOR_ALGORITHM(PanoramaIndoor, path);
    ROUTINE_FOR_ALGORITHM(YorkUrbanDB, path);
    ROUTINE_FOR_ALGORITHM(PanoContext, path);


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
        SaveToDisk("./cache/" + Tagify(path) + "_" + s + ".cereal", ts...);
    }

    template <class StringT, class ... Ts>
    inline void Load(const std::string & path, StringT && s, Ts & ... ts){
        LoadFromDisk("./cache/" + Tagify(path) + "_" + s + ".cereal", ts...);
    }




}

 
#endif