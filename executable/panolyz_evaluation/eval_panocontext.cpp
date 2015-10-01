#include "eval.hpp"

namespace panolyz {

    inline std::string PanoContextResultFilePath(const std::string & image) {
        return image + ".panocontext.mat";
    }

    struct ReconstructedModelPanoContext : ReconstructedModel {
        Image3ub image;
        std::vector<Polygon3> faces;
        std::vector<Vec3> vps;
        explicit ReconstructedModelPanoContext(const Image3ub & im, std::vector<Polygon3> && fs, std::vector<Vec3> && v) 
            : image(im), faces(std::move(fs)), vps(std::move(v)) {}
        virtual double depthAt(const Vec3 & direction) const {
            for (auto & f : faces) {
                auto inter = IntersectionOfLineAndPolygon(Ray3(Origin(), direction), f);
                if (inter.failed()) {
                    continue;
                }
                return Distance(inter.unwrap(), Origin());
            }
            std::cout << "Shape Not Closed At This Direction!" << std::endl;
            SHOULD_NEVER_BE_CALLED();
        }
        virtual void visualize(const std::vector<Vec3> & directions) const {
            gui::ResourceStore::set("texture", image);
            gui::SceneBuilder sb;
            sb.begin(faces).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();

            if (!directions.empty()) {
                sb.installingOptions().discretizeOptions.color = gui::Black;
                std::vector<Point3> testPoints;
                testPoints.reserve(directions.size());
                for (auto & dir : directions) {
                    double depth = depthAt(dir);
                    Point3 p = depth * normalize(dir);
                    testPoints.push_back(p);
                }
                sb.begin(testPoints).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).pointSize(10).end();
            }

           /* sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XPoints;
            sb.installingOptions().pointSize = 30;
            sb.add(gui::ColorAs(Point3(1, 0, 0), gui::Red));
            sb.add(gui::ColorAs(Point3(0, 1, 0), gui::Green));
            sb.add(gui::ColorAs(Point3(0, 0, 1), gui::Blue));
            sb.installingOptions().pointSize = 20;
            sb.add(gui::ColorAs(normalize(vps[0]) * .5, gui::Red));
            sb.add(gui::ColorAs(normalize(vps[1]) * .5, gui::Green));
            sb.add(gui::ColorAs(normalize(vps[2]) * .5, gui::Blue));
            sb.installingOptions().pointSize = 40;
            sb.add(gui::ColorAs(Origin(), gui::Black));*/
            sb.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(800, 800, Point2(400, 400), 400, Point3(1, 1, 1), Point3(0, 0, 0))));
        }
    };

    std::unique_ptr<ReconstructedModel> PredictionOfPanoContext(const std::string & impath, misc::Matlab & matlab) {

        Image3ub image = ImageRead(impath);
        misc::MAT panocontextResult(PanoContextResultFilePath(impath), misc::MAT::Read);
        auto result = panocontextResult.var("result");
        auto DET = result.field("DET");

        DenseMatd vpsData = result.field("vp");
        assert(vpsData.rows == 6 && vpsData.cols == 3);
        std::vector<Vec3> originalVps(3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                originalVps[i][j] = vpsData(i, j);
            }
        }

        std::vector<Vec3> vps = {
            { originalVps[2][1], originalVps[2][0], originalVps[2][2] },
            { originalVps[1][1], originalVps[1][0], originalVps[1][2] },
            { -originalVps[0][1], -originalVps[0][0], -originalVps[0][2] }
        };
        for (auto & vp : vps) {
            vp = normalize(vp);
        }
        auto trans = MakeMat4LocalToWorld(vps[0], vps[1], vps[2], Origin());

        int roomId = -1;
        for (int i = 0; i < DET.length(); i++) {
            std::string name = DET.field("name", i);
            if (name == "room") {
                roomId = i;
                break;
            }
        }
        if (roomId == -1) {
            std::cerr << "result.DET does not contain a room box!" << std::endl;
            return nullptr;
        }
        
        DenseMatd out_points_w = DET.field("out_points_w", roomId); // 8 x 3
        assert(out_points_w.cols == 3);

        std::vector<Point3> points(out_points_w.rows);
        double n = 0;
        for (int i = 0; i < out_points_w.rows; i++) {
            Vec4 c = trans * Vec4(out_points_w(i, 0), out_points_w(i, 1), out_points_w(i, 2), 1);
            points[i] = Vec3(c[0] / c[3], c[1] / c[3], c[2] / c[3]);
            n += points[i].ddot(points[i]);
        }
        for (auto & p : points) {
            p /= sqrt(n);
        }


        //   5 ---- 6
        //  /      /|
        // 8 ---- 7 |
        // | 4    | 3
        // |      |/
        // 1 ---- 2

        //   4 ---- 5
        //  /      /|
        // 7 ---- 6 |
        // | 3    | 2
        // |      |/
        // 0 ---- 1
        
        int faceCornerIds[6][4] = {
            {0, 1, 2, 3},
            {4, 5, 6, 7},
            {1, 0, 7, 6},
            {2, 1, 6, 5},
            {3, 2, 5, 4},
            {0, 3, 4, 7}
        };
        
        std::vector<Polygon3> faces;
        for (int i = 0; i < 6; i++) {
            Chain3 chain;
            for (int cid : faceCornerIds[i]) {
                chain.append(points[cid]);
            }
            faces.emplace_back(std::move(chain));
        }

        return std::make_unique<ReconstructedModelPanoContext>(image, std::move(faces), std::move(vps));

    }

}