#include "dataset.hpp"

namespace panoramix {
    namespace ml {        
        namespace annotations {

            void LoadFromDisk(const std::string & filename, sun::Panorama & panoInfo){
                std::ifstream in(filename);
                cereal::JSONInputArchive archive(in);
                archive(cereal::make_nvp("cameras", panoInfo.cameras),
                    cereal::make_nvp("objects", panoInfo.objects),
                    cereal::make_nvp("PointCloud", panoInfo.pointCloud));
                std::cout << "file \"" << filename << "\" loaded" << std::endl;
                in.close();
            }

            void SaveToDisk(const std::string & filename, const sun::Panorama & panoInfo){
                std::ofstream out(filename);
                cereal::JSONOutputArchive archive(out);
                archive(cereal::make_nvp("cameras", panoInfo.cameras),
                    cereal::make_nvp("objects", panoInfo.objects),
                    cereal::make_nvp("PointCloud", panoInfo.pointCloud));
                std::cout << "file \"" << filename << "\" saved" << std::endl;
                out.close();
            }

        }

    }
}
 