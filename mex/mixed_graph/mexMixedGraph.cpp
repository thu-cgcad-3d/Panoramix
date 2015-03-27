#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/gui/visualize2d.hpp"

using namespace panoramix;
using namespace panoramix::core;

// returns (fitted plane, inlier ratio)
std::pair<Plane3, double> FitPlane(const std::vector<Point3> & pcloud){
    //DenseMatd mean;
    //cv::PCA pca(pcloud, mean, CV_PCA_DATA_AS_COL);

    //double pToPlaneThres = pca.eigenvalues.at<double>(2);
    //int numInliners = 0;

    //Vec3 norm = normalize(pca.eigenvectors.row(2));
    //Point3 x0 = pca.mean;


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

struct MixedGraphGroundTruthTable {
    Imagef depths;
    HandledTable<RegionHandle, RegionGroundTruthData> regionGT;
    HandledTable<RegionBoundaryHandle, RegionBoundaryGroundTruthData> boundaryGT;

    RegionGroundTruthData & operator[](RegionHandle rh) { return regionGT[rh]; }
    const RegionGroundTruthData & operator[](RegionHandle rh) const { return regionGT[rh]; }

    RegionBoundaryGroundTruthData & operator[](RegionBoundaryHandle bh) { return boundaryGT[bh]; }
    const RegionBoundaryGroundTruthData & operator[](RegionBoundaryHandle bh) const { return boundaryGT[bh]; }

    explicit MixedGraphGroundTruthTable(const MixedGraph & mg) 
        : regionGT(mg.internalComponents<RegionData>().size()), 
        boundaryGT(mg.internalConstraints<RegionBoundaryData>().size()) {}
};

struct MG {
    core::View<core::PerspectiveCamera> view;
    std::vector<core::Classified<core::Line2>> lines;
    std::vector<core::Vec3> vps;
    core::Imagei segmentedImage;
    std::vector<RegionHandle> regionHandles;

    MixedGraph g;
    MixedGraphPropertyTable p;

    std::shared_ptr<MixedGraphGroundTruthTable> gt;

    ~MG() { printf("BoundaryGraph destroyed\n"); }
};


MG CreateMG(const Image & image, double f, const Point2 & c) {

    core::View<core::PerspectiveCamera> view;
    std::vector<core::Classified<core::Line2>> lines;
    std::vector<core::Vec3> vps;

    core::MixedGraph mg;

    core::Imagei segmentedImage;
    core::MixedGraphPropertyTable props;

    core::LineSegmentExtractor lineExtractor;
    lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
    lines = core::ClassifyEachAs(lineExtractor(image, 2), -1);

    view.camera = core::PerspectiveCamera(image.cols, image.rows, c, f, 
        core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1));
    view.image = image;

    // classify lines
    vps = core::EstimateVanishingPointsAndClassifyLines(view.camera, lines);
    // remove non-manhattan vps
    vps.erase(vps.begin() + 3, vps.end());
    for (auto & l : lines){
        if (l.claz >= 3){
            l.claz = -1;
        }
    }

    // append lines
    core::AppendLines(mg, lines, view.camera, vps);

    // append regions
    core::SegmentationExtractor segmenter;
    segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
    segmenter.params().c = 100.0;
    std::vector<core::Line2> pureLines(lines.size());
    for (int i = 0; i < lines.size(); i++)
        pureLines[i] = lines[i].component;
    int segmentsNum = 0;
    std::tie(segmentedImage, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

    auto regionHandles = core::AppendRegions(mg, segmentedImage, view.camera, 0.001, 0.001, 3, 1);

    // attach constraints
    props = core::MakeMixedGraphPropertyTable(mg, vps);
    core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 120.0);
    core::AttachWallConstriants(mg, props, M_PI / 100.0);
    core::ResetVariables(mg, props);

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
    std::vector<gui::Color> colors(g.regionHandles.size(), gui::ColorTag::Yellow);
    std::vector<gui::Color> oColors = { gui::ColorTag::Red, gui::ColorTag::Green, gui::ColorTag::Blue };
    std::vector<gui::Color> noColors = { gui::ColorTag::DimGray, gui::ColorTag::Gray, gui::ColorTag::DarkGray };
    for (int i = 0; i < g.regionHandles.size(); i++){
        if (g.regionHandles[i].invalid())
            continue;
        auto rh = g.regionHandles[i];
        auto & p = props[rh];
        auto & color = colors[i];
        if (!p.used){
            color = gui::ColorTag::Black;
        }
        else{
            if (p.orientationClaz != -1){
                color = oColors[p.orientationClaz];
            }
            else if (p.orientationNotClaz != -1){
                color = noColors[p.orientationNotClaz];
            }
            else{
                color = gui::ColorTag::White;
            }
        }
    }
    gui::Visualizer2D::Params params;
    params.colorTable = gui::ColorTable(colors, gui::ColorTag::White);
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
    g.gt = std::make_shared<MixedGraphGroundTruthTable>(g.g);
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
                cv::line(im, p1, p2, gui::Color(gt[b.topo.hd].ot == R1OccludesR2 ? gui::ColorTag::Blue : gui::ColorTag::Red), 1);
            }
        }
    }
    cv::imshow("GT Occlusions", im);
    cv::waitKey();
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
        printf("MixedGraph created\n");
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

    if (cmd == "showOc"){
        if (!g.gt){
            mexErrMsgTxt("No GroundTruth depths installed yet");
            return;
        }

        ShowGTOcclusions(g);
        return;
    }


    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
