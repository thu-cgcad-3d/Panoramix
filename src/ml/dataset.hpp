#ifndef PANORAMIX_ML_DATASET_HPP
#define PANORAMIX_ML_DATASET_HPP

#include "../core/basic_types.hpp"
 
namespace panoramix {
    namespace ml {

        // dataset annotations
        namespace annotations {
            
            // SUN dataset
            namespace sun {

                struct Rotation {
                    double pitch, yaw, roll;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(pitch), CEREAL_NVP(yaw), CEREAL_NVP(roll));
                    }
                };

                struct Point2 {
                    double x, y;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(x), CEREAL_NVP(y));
                    }
                };

                struct Point3 {
                    template <class Archive>
                    void serialize(Archive & ar) {}
                };

                struct Size {
                    double width, height;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(width), CEREAL_NVP(height));
                    }
                };

                struct Image {
                    std::string folder, file;
                    int frame;
                    Point2 principle;
                    Size resolution;
                    Size film;
                    std::string type;
                    std::string plane;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(folder), CEREAL_NVP(file), CEREAL_NVP(frame),
                            CEREAL_NVP(principle), CEREAL_NVP(resolution), CEREAL_NVP(film),
                            CEREAL_NVP(type), CEREAL_NVP(plane));
                    }
                };

                struct Camera {
                    std::vector<double> modelview;
                    Rotation rotation;
                    double focal_length;
                    std::vector<Image> images;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(modelview),
                            CEREAL_NVP(rotation),
                            CEREAL_NVP(focal_length),
                            CEREAL_NVP(images));
                    }
                };

                struct Projection {
                    int camera;
                    int image;
                    Point2 position2D;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(camera),
                            CEREAL_NVP(image),
                            CEREAL_NVP(position2D));
                    }
                };

                struct ObjectPoint{
                    Point3 position3D;
                    std::vector<Projection> projection;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(position3D),
                            CEREAL_NVP(projection));
                    }
                };

                struct Object {
                    int type;
                    std::string name;
                    std::string time;
                    std::string creator;
                    std::vector<ObjectPoint> points;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(type), CEREAL_NVP(name),
                            CEREAL_NVP(time), CEREAL_NVP(creator),
                            CEREAL_NVP(points));
                    }
                };

                struct PointCloud {
                    int number;
                    std::string vertex;
                    std::string color;
                    template <class Archive>
                    void serialize(Archive & ar) {
                        ar(CEREAL_NVP(number),
                            CEREAL_NVP(vertex),
                            CEREAL_NVP(color));
                    }
                };


                // the panorama annotation
                struct Panorama {
                    std::vector<Camera> cameras;
                    std::vector<Object> objects;
                    PointCloud pointCloud;
                };
            }


            void LoadFromDisk(const std::string & filename, sun::Panorama & panoInfo);
            void SaveToDisk(const std::string & filename, const sun::Panorama & panoInfo);

        }

    }
}
 
 
#endif