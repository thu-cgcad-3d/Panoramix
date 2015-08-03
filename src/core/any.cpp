#include <unordered_map>
#include "any.hpp"

namespace pano {
    namespace core {

        static std::unordered_map<std::string, Any> staticDict;


        void StaticStorage::set(const std::string & name, const Any & val) {
            staticDict[name] = val;
        }

        void StaticStorage::set(const std::string & name, Any && val) {
            staticDict[name] = std::move(val);
        }


        bool StaticStorage::has(const std::string & name) {
            return staticDict.find(name) != staticDict.end();
        }

        Any StaticStorage::get(const std::string & name) {
            return staticDict.at(name);
        }

    }
}