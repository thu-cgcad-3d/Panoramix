#pragma once

#include "../core/serialization.hpp"
#include "../core/utility.hpp"
#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        //// annotation
        //struct AnnotedPolygon {
        //    Polygon3 polygon;
        //    SegControl control;
        //    template <class Archiver>
        //    inline void serialize(Archiver & ar) {
        //        ar(polygon, control);
        //    }
        //};
        //struct AnnotedOcclusion {
        //    Chain3 chain; // left is always front!
        //    template <class Archiver>
        //    inline void serialize(Archiver & ar) {
        //        ar(chain);
        //    }
        //};

        //struct PIAnnotation {
        //    Image originalImage;
        //    Image rectifiedImage;
        //    bool extendedOnTop;
        //    bool extendedOnBottom;

        //    PanoramicView view;
        //    std::vector<Vec3> vps;
        //    int vertVPId;
        //    std::vector<Classified<Line3>> lines;
        //    std::vector<AnnotedPolygon> polygons;
        //    std::vector<AnnotedOcclusion> occlusions;
        //    std::vector<Line3> foldingLines;

        //    PIAnnotation() : vertVPId(-1) {}

        //    template <class Archiver>
        //    inline void load(Archiver & ar, std::int32_t version) {
        //        if (version == 0) {
        //            ar(originalImage, rectifiedImage, view, vps, vertVPId, lines, polygons, occlusions);
        //        } else if(version == 1){
        //            ar(originalImage, rectifiedImage, view, vps, vertVPId, lines, polygons, occlusions, foldingLines);
        //        } else {
        //            NOT_IMPLEMENTED_YET();
        //        }
        //    }

        //    template <class Archiver>
        //    inline void save(Archiver & ar, std::int32_t version) const {
        //        if (version == 0) {
        //            ar(originalImage, rectifiedImage, view, vps, vertVPId, lines, polygons, occlusions);
        //        } else if (version == 1) {
        //            ar(originalImage, rectifiedImage, view, vps, vertVPId, lines, polygons, occlusions, foldingLines);
        //        } else {
        //            NOT_IMPLEMENTED_YET();
        //        }
        //    }
        //};


        //std::string AnnotationFilePath(const std::string & imagePath);

        //PIAnnotation LoadOrInitializeNewAnnotation(const std::string & imagePath);

        //void EditAnnotation(PIAnnotation & anno);

        //void SaveAnnotation(const std::string & imagePath, const PIAnnotation & anno);

        //// AttachAnnotatedPolygonsAndOcclusions
        //void AttachAnnotatedPolygonsAndOcclusions(PIGraph & mg,
        //    const std::vector<AnnotedPolygon> & polygons, const std::vector<AnnotedOcclusion> & occs,
        //    double polygonBoundarySampleStepAngle = DegreesToRadians(1),
        //    double occChainSampleStepAngle = DegreesToRadians(0.1),
        //    double occChainToBndPieceAngleThres = DegreesToRadians(3));



        

        struct PILayoutAnnotation {
            Image originalImage;
            Image rectifiedImage;
            bool extendedOnTop;
            bool extendedOnBottom;

            PanoramicView view;
            std::vector<Vec3> vps;
            int vertVPId;
                        
            std::vector<Vec3> corners;
            std::map<std::pair<int, int>, int> corners2border;
            std::vector<std::pair<int, int>> border2corners; // from, to
            std::vector<bool> border2connected;
            std::vector<std::pair<int, int>> border2faces; // left, right

            int addBorder(int c1, int c2);
            int splitBorderBy(int b, int c);

            std::vector<std::vector<int>> face2corners;
            std::vector<SegControl> face2control;
            std::vector<Plane3> face2plane;
            std::vector<std::pair<int, int>> coplanarFaces;
            std::vector<Polygon3> clutters;

            PILayoutAnnotation() : vertVPId(-1) {}
            int ncorners() const { return corners.size(); }
            int nborders() const { return border2corners.size(); }
            int nfaces() const { return face2corners.size(); }
            void generateFacesWithBorders();

            template <class Archiver>
            inline void serialize(Archiver & ar, std::int32_t version) {
                if (version == 0) {
                    ar(originalImage, rectifiedImage, extendedOnTop, extendedOnBottom);
                    ar(view, vps, vertVPId);
                    ar(corners, corners2border, border2corners, border2connected, border2faces);
                    ar(face2corners, face2control, face2plane, coplanarFaces, clutters);
                } else {
                    NOT_IMPLEMENTED_YET();
                }
            }
        };

        std::string LayoutAnnotationFilePath(const std::string & imagePath);

        PILayoutAnnotation LoadOrInitializeNewLayoutAnnotation(const std::string & imagePath);

        void EditLayoutAnnotation(PILayoutAnnotation & anno);

        void SaveLayoutAnnotation(const std::string & imagePath, const PILayoutAnnotation & anno);

    }
}


//CEREAL_CLASS_VERSION(pano::experimental::PIAnnotation, 0);

CEREAL_CLASS_VERSION(pano::experimental::PILayoutAnnotation, 0);