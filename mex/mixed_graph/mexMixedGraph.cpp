#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/gui/visualize2d.hpp"

using namespace panoramix;
using namespace panoramix::core;


struct MG {
    core::View<core::PerspectiveCamera> view;
    std::vector<core::Classified<core::Line2>> lines;
    std::vector<core::Vec3> vps;
    core::Imagei segmentedImage;

    MixedGraph g;
    MixedGraphPropertyTable p;
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

    gui::Visualizer2D::Params vParams;
    vParams.colorTable = gui::CreateRandomColorTableWithSize(segmentsNum);
    gui::Visualizer2D(segmentedImage, vParams)
        << gui::manip2d::Show();

    core::AppendRegions(mg, segmentedImage, view.camera, 0.001, 0.001, 3, 1);

    // attach constraints
    props = core::MakeMixedGraphPropertyTable(mg, vps);
    core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 120.0);
    core::AttachWallConstriants(mg, props, M_PI / 100.0);
    core::ResetVariables(mg, props);

    return MG{ std::move(view), std::move(lines), std::move(vps), std::move(segmentedImage), std::move(mg), std::move(props) };

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
    std::vector<gui::Color> colors(mg.internalComponents<core::RegionData>().size());
    std::vector<gui::Color> oColors = { gui::ColorTag::Red, gui::ColorTag::Green, gui::ColorTag::Blue };
    std::vector<gui::Color> noColors = { gui::ColorTag::DimGray, gui::ColorTag::Gray, gui::ColorTag::DarkGray };
    for (auto & r : mg.components<core::RegionData>()){
        auto & p = props[r.topo.hd];
        auto & color = colors[r.topo.hd.id];
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Get the command string
    char cmdbuffer[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmdbuffer, sizeof(cmdbuffer)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    std::string cmd = cmdbuffer;

    const mxArray ** argv = prhs + 1;
    int argc = nrhs - 1;

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
        return;
    }
    
    if (cmd == "delete") {
        destroyObject<MG>(prhs[1]);
        printf("BoundaryGraph destroyed\n");
        return;
    }
    
    // Get the class instance pointer from the second input
    MG & g = *convertMat2Ptr<MG>(prhs[1]);

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


    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
