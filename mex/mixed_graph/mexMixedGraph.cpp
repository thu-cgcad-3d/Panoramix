#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/algorithms.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

using namespace panoramix;
using namespace panoramix::core;
using namespace panoramix::experimental;

namespace panoramix {
    namespace core{
        using namespace experimental;
    }
}

// returns (fitted plane, inlier ratio)
std::pair<Plane3, double> FitPlane(const std::vector<Point3> & pcloud){

    auto results = EigenVectorAndValuesFromPoints(pcloud);
    std::sort(results.begin(), results.end(), std::greater<>());
    double pToPlaneThres = results[2].score;
    int numInliners = 0;

    Vec3 norm = normalize(results[2].component);
    Point3 x0(0, 0, 0);
    for (auto & p : pcloud){
        x0 += p;
    }
    x0 /= double(pcloud.size());

    for (int i = 0; i < pcloud.size(); i++) {
        Vec3 w = pcloud[i] - x0;
        double D = abs(norm.dot(w));
        if (D < pToPlaneThres)
            numInliners ++;
    }
    return std::make_pair(Plane3(x0, norm), numInliners / double(pcloud.size()));
}

inline double DepthAt(const Vec3 & direction, const Plane3 & plane){
    Ray3 ray(Point3(0, 0, 0), direction);
    return norm(IntersectionOfLineAndPlane(ray, plane).position);
}

inline double DepthAt(const Vec3 & direction, const Line3 & line){
    Ray3 ray(Point3(0, 0, 0), direction);
    return norm(DistanceBetweenTwoLines(ray, line.infiniteLine()).second.first);
}



enum OcclusionType {
    R1OccludesR2 = 0,
    R2OccludesR1,
    Connected
};

struct RegionGroundTruthData {
    Plane3 fittedPlane;
    double inlinerRatio;
};

struct RegionBoundaryGroundTruthData {
    OcclusionType ot;
};

struct RLGraphGroundTruthTable {
    Imagef depths;
    HandledTable<RegionHandle, RegionGroundTruthData> regionGT;
    HandledTable<RegionBoundaryHandle, RegionBoundaryGroundTruthData> boundaryGT;

    RegionGroundTruthData & operator[](RegionHandle rh) { return regionGT[rh]; }
    const RegionGroundTruthData & operator[](RegionHandle rh) const { return regionGT[rh]; }

    RegionBoundaryGroundTruthData & operator[](RegionBoundaryHandle bh) { return boundaryGT[bh]; }
    const RegionBoundaryGroundTruthData & operator[](RegionBoundaryHandle bh) const { return boundaryGT[bh]; }

    explicit RLGraphGroundTruthTable(const RLGraph & mg) 
        : regionGT(mg.internalComponents<RegionData>().size()), 
        boundaryGT(mg.internalConstraints<RegionBoundaryData>().size()) {}
};

struct MG {
    core::View<core::PerspectiveCamera> view;
    std::vector<core::Classified<core::Line2>> lines;
    std::vector<core::Vec3> vps;
    core::Imagei segmentedImage;
    std::vector<RegionHandle> regionHandles;

    RLGraph g;
    RLGraphPropertyTable p;

    std::shared_ptr<RLGraphGroundTruthTable> gt;

    ~MG() { printf("BoundaryGraph destroyed\n"); }
};


MG CreateMG(const Image & image, double f, const Point2 & c) {

    core::View<core::PerspectiveCamera> view;
    std::vector<core::Classified<core::Line2>> lines;
    std::vector<core::Vec3> vps;

    RLGraph mg;

    core::Imagei segmentedImage;
    RLGraphPropertyTable props;

    core::LineSegmentExtractor lineExtractor;
    lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
    lines = core::ClassifyEachAs(lineExtractor(image, 2), -1);

    view.camera = core::PerspectiveCamera(image.cols, image.rows, c, f, 
        core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1));
    view.image = image;

    // classify lines
    vps = EstimateVanishingPointsAndClassifyLines(view.camera, lines);
    // remove non-manhattan vps
    vps.erase(vps.begin() + 3, vps.end());
    for (auto & l : lines){
        if (l.claz >= 3){
            l.claz = -1;
        }
    }

    // append lines
    AppendLines(mg, lines, view.camera, vps);

    // append regions
    core::SegmentationExtractor segmenter;
    segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
    segmenter.params().c = 100.0;
    std::vector<core::Line2> pureLines(lines.size());
    for (int i = 0; i < lines.size(); i++)
        pureLines[i] = lines[i].component;
    int segmentsNum = 0;
    std::tie(segmentedImage, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

    auto regionHandles = AppendRegions(mg, segmentedImage, view.camera, 0.001, 0.001, 3, 1);

    // attach constraints
    props = MakeRLGraphPropertyTable(mg, vps);
    AttachPrincipleDirectionConstraints(mg, props, M_PI / 120.0);
    AttachWallConstriants(mg, props, M_PI / 100.0);
    ResetVariables(mg, props);

    return MG{ std::move(view), std::move(lines), std::move(vps), 
        std::move(segmentedImage), 
        std::move(regionHandles), 
        std::move(mg), std::move(props) };

}


void ShowSegmentations(const MG & g){
    int segmentsNum = MinMaxValOfImage(g.segmentedImage).second + 1;
    gui::Visualizer2D::Params vParams;
    vParams.colorTable = gui::CreateRandomColorTableWithSize(segmentsNum);
    gui::Visualizer2D(g.segmentedImage, vParams)
        << gui::manip2d::Show();
}


void ShowVPsAndLines(const MG & g){
    std::vector<core::Classified<core::Ray2>> vpRays;
    for (int i = 0; i < 3; i++){
        std::cout << "vp[" << i << "] = " << g.view.camera.screenProjection(g.vps[i]) << std::endl;
        for (double a = 0; a <= M_PI * 2.0; a += 0.1){
            core::Point2 p = core::Point2(g.view.image.cols / 2, g.view.image.rows / 2) + core::Vec2(cos(a), sin(a)) * 1000.0;
            vpRays.push_back(core::ClassifyAs(core::Ray2(p, (g.view.camera.screenProjectionInHPoint(g.vps[i]) - core::HPoint2(p, 1.0)).numerator), i));
        }
    }
    gui::Visualizer2D(g.view.image)
        << gui::manip2d::SetColorTable(gui::ColorTable(gui::ColorTableDescriptor::RGB).appendRandomizedGreyColors(g.vps.size() - 3))
        << gui::manip2d::SetThickness(1)
        << vpRays
        << gui::manip2d::SetThickness(2)
        << g.lines
        << gui::manip2d::Show(0);
}

void ShowRegionOrientationConstraints(const MG & g){
    auto & mg = g.g;
    auto & props = g.p;
    std::vector<gui::Color> colors(g.regionHandles.size(), gui::Yellow);
    std::vector<gui::Color> oColors = { gui::Red, gui::Green, gui::Blue };
    std::vector<gui::Color> noColors = { gui::DimGray, gui::Gray, gui::DarkGray };
    for (int i = 0; i < g.regionHandles.size(); i++){
        if (g.regionHandles[i].invalid())
            continue;
        auto rh = g.regionHandles[i];
        auto & p = props[rh];
        auto & color = colors[i];
        if (!p.used){
            color = gui::Black;
        }
        else{
            if (p.orientationClaz != -1){
                color = oColors[p.orientationClaz];
            }
            else if (p.orientationNotClaz != -1){
                color = noColors[p.orientationNotClaz];
            }
            else{
                color = gui::White;
            }
        }
    }
    gui::Visualizer2D::Params params;
    params.colorTable = gui::ColorTable(colors, gui::White);
    gui::Visualizer2D(g.segmentedImage, params)
        << gui::manip2d::Show();
}

Imagef ComputeDepths(const MG & g){
    Imagef depths(g.view.image.size(), 0.0);
    for (auto it = depths.begin(); it != depths.end(); ++it){
        auto p = it.pos();
        Vec3 dir = g.view.camera.spatialDirection(p);
        int regionId = g.segmentedImage(p);
        Plane3 plane = Instance(g.g, g.p, RegionHandle(regionId));
        double d = norm(IntersectionOfLineAndPlane(Ray3(g.view.camera.eye(), dir), plane).position);
        *it = d;
    }
    return depths;
}


void InstallGT(MG & g, const Imagef & depths){
    g.gt = std::make_shared<RLGraphGroundTruthTable>(g.g);
    auto & gt = *g.gt;
    gt.depths = depths;
    auto & mg = g.g;
    auto & segs = g.segmentedImage;

    double maxDepths, minDepths;
    std::tie(minDepths, maxDepths) = MinMaxValOfImage(depths);

    auto regionPoints = mg.createComponentTable<RegionData, std::vector<Point3>>();
    for (auto it = segs.begin(); it != segs.end(); ++it){
        double d = depths(it.pos());
        RegionHandle rh = g.regionHandles[*it];
        if (rh.invalid())
            continue;
        regionPoints[rh].push_back(normalize(g.view.camera.spatialDirection(it.pos())) * d);
    }

    // fit planes
    for (auto & r : mg.components<core::RegionData>()){
        auto & rgt = gt[r.topo.hd];
        std::tie(rgt.fittedPlane, rgt.inlinerRatio) = FitPlane(regionPoints[r.topo.hd]);
    }
    // recognize occlusions
    for (auto & b : mg.constraints<core::RegionBoundaryData>()){
        auto rh1 = b.topo.component<0>();
        auto rh2 = b.topo.component<1>();
        auto & plane1 = gt[rh1].fittedPlane;
        auto & plane2 = gt[rh2].fittedPlane;
        double depthsSum1 = 0.0, depthsSum2 = 0.0;
        int num = 0;
        for (auto & cs : b.data.normalizedEdges){
            for (auto & c : cs){
                double d1 = gt[rh1].inlinerRatio < 0.5 ?
                    depths(ToPixelLoc(g.view.camera.screenProjection(mg.data(rh1).normalizedCenter))) :
                    norm(IntersectionOfLineAndPlane(Ray3(g.view.camera.eye(), c), plane1).position);
                depthsSum1 += d1;
                double d2 = gt[rh2].inlinerRatio < 0.5 ?
                    depths(ToPixelLoc(g.view.camera.screenProjection(mg.data(rh2).normalizedCenter))) :
                    norm(IntersectionOfLineAndPlane(Ray3(g.view.camera.eye(), c), plane2).position);
                depthsSum2 += d2;
                num++;
            }
        }
        double dmean1 = depthsSum1 / num, dmean2 = depthsSum2 / num;
        double thres = (depthsSum1 + depthsSum2) / num / 50.0;
        if (dmean1 - dmean2 > thres){
            gt[b.topo.hd].ot = R2OccludesR1;
        }
        else if (dmean2 - dmean1 > thres){
            gt[b.topo.hd].ot = R1OccludesR2;
        }
        else{
            gt[b.topo.hd].ot = Connected;
        }
    }
}


void ShowGTOcclusions(const MG & g){
    auto & mg = g.g;
    auto & props = g.p;
    auto & gt = *g.gt;

    // render disconnected boundaries
    Image im = g.view.image.clone();
    for (auto & b : mg.constraints<core::RegionBoundaryData>()){
        if (gt[b.topo.hd].ot == Connected)
            continue;
        for (auto & e : b.data.normalizedEdges){
            if (e.size() <= 1) continue;
            for (int i = 1; i < e.size(); i++){
                auto p1 = core::ToPixelLoc(g.view.camera.screenProjection(e[i - 1]));
                auto p2 = core::ToPixelLoc(g.view.camera.screenProjection(e[i]));
                cv::clipLine(cv::Rect(0, 0, im.cols, im.rows), p1, p2);
                cv::line(im, p1, p2, gui::Color(gt[b.topo.hd].ot == R1OccludesR2 ? gui::Blue : gui::Red), 1);
            }
        }
    }
    cv::imshow("GT Occlusions", im);
    cv::waitKey();
}

void ShowGTPlanes(const MG & g){
    auto & mg = g.g;
    auto & gt = *g.gt;

    struct ComponentID {
        int handleID;
        bool isRegion;
    };
    gui::ResourceStore::set("texture", g.view.sampled(PanoramicCamera(250, Point3(0, 0, 0), Point3(1, 0, 0), Vec3(0, 0, 1))).image);
    gui::Visualizer viz("mixed graph optimizable");
    viz.renderOptions.bwColor = 1.0;
    viz.renderOptions.bwTexColor = 0.0;
    viz.installingOptions.discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
    std::vector<std::pair<ComponentID, gui::Colored<gui::SpatialProjectedPolygon>>> spps;
    std::vector<gui::Colored<core::Line3>> lines;

    for (auto & c : mg.components<RegionData>()){
        auto uh = c.topo.hd;
        auto & region = c.data;
        gui::SpatialProjectedPolygon spp;
        // filter corners
        core::ForeachCompatibleWithLastElement(c.data.normalizedContours.front().begin(), c.data.normalizedContours.front().end(),
            std::back_inserter(spp.corners),
            [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
            return core::AngleBetweenDirections(a, b) > M_PI / 1000.0;
        });
        if (spp.corners.size() < 3)
            continue;

        spp.projectionCenter = core::Point3(0, 0, 0);
        spp.plane = gt[uh].fittedPlane;
        if (spp.plane.distanceTo(Point3(0, 0, 0)) < 0.2)
            continue;
        assert(!HasValue(spp.plane, IsInfOrNaN<double>));
        spps.emplace_back(ComponentID{ uh.id, true }, std::move(gui::ColorAs(spp, gui::Transparent)));
    }

    viz.begin(spps).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
    viz.installingOptions.discretizeOptions.color = gui::DarkGray;
    viz.installingOptions.lineWidth = 5.0;
    viz.add(lines);

   /* std::vector<core::Line3> connectionLines;
    for (auto & c : mg.constraints<RegionBoundaryData>()){
        auto & samples = c.data.normalizedSampledPoints;
        auto inst1 = gt[c.topo.component<0>()].fittedPlane;
        auto inst2 = gt[c.topo.component<1>()].fittedPlane;
        for (auto & ss : samples){
            for (auto & s : ss){
                double d1 = DepthAt(s, inst1);
                double d2 = DepthAt(s, inst2);
                connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
            }
        }
    }*/

    viz.installingOptions.discretizeOptions.color = gui::Black;
    viz.installingOptions.lineWidth = 1.0;
    //viz.add(connectionLines);

    viz.renderOptions.renderMode = gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines;
    viz.renderOptions.backgroundColor = gui::White;
    viz.renderOptions.bwColor = 0.5;
    viz.renderOptions.bwTexColor = 0.5;
    viz.camera(core::PerspectiveCamera(1000, 800, core::Point2(500, 400), 800, Point3(-1, 1, 1), Point3(0, 0, 0)));
    viz.show(true, false);

    gui::ResourceStore::clear();
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Get the command string
    char cmdbuffer[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmdbuffer, sizeof(cmdbuffer)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    std::string cmd = cmdbuffer;

    const mxArray ** argv = prhs + 1;
    int argc = nrhs - 1;
    int outc = nlhs;
    mxArray ** outv = plhs;

    if (cmd == "new") {
        if (argc != 3){
            mexErrMsgTxt("Wrong Arguments Num.");
            return;
        }
        Image image;
        misc::Matlab::GetVariable(argv[0], image);
        double focal = mxGetScalar(argv[1]);
        Point2 pp;
        misc::Matlab::GetVariable(argv[2], pp);

        if (image.channels() == 3)
            cv::cvtColor(image, image, CV_BGR2RGB);

        plhs[0] = convertPtr2Mat(new MG(CreateMG(image, focal, pp)));
        printf("RLGraph created\n");
        return;
    }
    
    if (cmd == "delete") {
        destroyObject<MG>(prhs[1]);
        return;
    }
    
    // Get the class instance pointer from the second input
    MG & g = *convertMat2Ptr<MG>(prhs[1]);
    argc--;
    argv++;

    if (cmd == "showSegs"){
        ShowSegmentations(g);
        return;
    }

    if (cmd == "showVPsLines"){
        ShowVPsAndLines(g);
        return;
    }

    if (cmd == "showROC"){
        ShowRegionOrientationConstraints(g);
        return;
    }

    if (cmd == "solve"){
        auto & mg = g.g;
        auto & props = g.p;
        core::SolveVariablesUsingInversedDepths(mg, props, false);
        core::NormalizeVariables(mg, props);
        std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
        return;
    }

    if (cmd == "looseC"){
        core::LooseOrientationConstraintsOnComponents(g.g, g.p, 0.2, 0, 0.05);
        return;
    }

    if (cmd == "show"){
        core::Visualize(g.view, g.g, g.p);
        return;
    }

    if (cmd == "depths"){
        if (outc == 0){
            return;
        }
        if (outc > 1){
            mexErrMsgTxt("Too Many Outputs");
            return;
        }

        Imagef depths = ComputeDepths(g);
        mxArray * depthsMXA = static_cast<mxArray*>(misc::Matlab::PutVariable(depths));
        outv[0] = depthsMXA;
        return;
    }


    if (cmd == "installGT"){
        if (argc != 1){
            mexErrMsgTxt("Wrong Arguments Num");
            return;
        }

        Imagef depths;
        misc::Matlab::GetVariable(argv[0], depths);
        InstallGT(g, depths);
        return;
    }

    if (cmd == "depthsGT"){
        if (outc == 0){
            return;
        }
        if (outc > 1){
            mexErrMsgTxt("Too Many Outputs");
            return;
        }
        if (!g.gt){
            mexErrMsgTxt("No GroundTruth depths installed yet");
            return;
        }
        mxArray * depthsMXA = static_cast<mxArray*>(misc::Matlab::PutVariable(g.gt->depths));
        outv[0] = depthsMXA;
        return;
    }



    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
