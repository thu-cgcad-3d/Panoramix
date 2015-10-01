#include "eval.hpp"

namespace panolyz {
    
    struct GTModel : ReconstructedModel{
        PILayoutAnnotation anno;
        std::vector<Polygon3> polygons;
        explicit GTModel(PILayoutAnnotation && a) : anno(std::move(a)), polygons(anno.nfaces()) {
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & plane = anno.face2plane[i];
                auto & poly = polygons[i];
                poly.normal = plane.normal;
                for (int c : anno.face2corners[i]) {
                    Ray3 ray(Origin(), anno.corners[c]);
                    poly.corners.push_back(IntersectionOfLineAndPlane(ray, plane).position);
                }
            }
        }
        virtual double depthAt(const Vec3 & direction) const override {
            Ray3 ray(Origin(), direction);
            double depth = std::numeric_limits<double>::infinity();
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & poly = polygons[i];
                auto inter = IntersectionOfLineAndPolygon(ray, poly);
                if (inter.failed()) {
                    continue;
                }
                double curDepth = norm(inter.unwrap());
                if (curDepth < depth) {
                    depth = curDepth;
                }
            }
            return depth;
        }
        virtual bool isValid(const Vec3 & direction) const override {
            Ray3 ray(Origin(), direction);
            for (auto & clutterPoly : anno.clutters) {
                auto inter = IntersectionOfLineAndPolygon(ray, clutterPoly);
                if (!inter.failed()) {
                    return false;
                }
            }
            return true;
        }
        virtual void visualize(const std::vector<Vec3> & directions) const override {
            gui::SceneBuilder viz;
            viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
            std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int>> spps;

            gui::ResourceStore::set("texture", anno.rectifiedImage);
            for (int face = 0; face < anno.nfaces(); face++) {

                gui::SpatialProjectedPolygon spp;
                spp.projectionCenter = core::Point3(0, 0, 0);
                spp.plane = anno.face2plane[face];
                for (int c : anno.face2corners[face]) {
                    spp.corners.push_back(normalize(anno.corners[c]));
                }

                static const gui::ColorTable rgbTable = gui::RGB;
                auto & control = anno.face2control[face];

                spps.push_back(core::DecorateAs(std::move(gui::ColorAs(spp,
                    control.dof() == 1
                    ? rgbTable[control.orientationClaz]
                    : (control.dof() == 2
                    ? rgbTable[control.orientationNotClaz]
                    : gui::Color(gui::White)))), face));
            }

            viz.begin(spps/*, sppCallbackFun*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();

            if (!directions.empty()) {
                viz.installingOptions().discretizeOptions.color = gui::Black;
                std::vector<Point3> testPoints;
                testPoints.reserve(directions.size());
                for (auto & dir : directions) {
                    double depth = depthAt(dir);
                    Point3 p = depth * normalize(dir);
                    testPoints.push_back(p);
                }
                viz.begin(testPoints).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).pointSize(10).end();
            }

            viz.show(true, false, gui::RenderOptions().cullFrontFace(true).cullBackFace(false).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
        }
    };

    std::unique_ptr<ReconstructedModel> PredictionOfGT(const std::string & impath) {
        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
        return std::make_unique<GTModel>(std::move(anno));
    }

}

