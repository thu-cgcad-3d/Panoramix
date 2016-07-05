#include <chrono>

#include "segmentation.hpp"
#include "geo_context.hpp"
#include "line_detection.hpp"
#include "panoramix.hpp"

template <class T> double ElapsedInMS(const T &start) {
  return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
             std::chrono::system_clock::now() - start)
      .count();
}

const std::string PanoramixOptions::parseOption(bool b) {
  return b ? "_on" : "_off";
}

std::string PanoramixOptions::algorithmOptionsTag() const {
  std::stringstream ss;
  ss << "_LayoutVersion" << LayoutVersion << parseOption(useWallPrior)
     << parseOption(usePrincipleDirectionPrior)
     << parseOption(useGeometricContextPrior) << parseOption(useGTOcclusions)
     << parseOption(looseLinesSecondTime) << parseOption(looseSegsSecondTime)
     << parseOption(restrictSegsSecondTime);
  if (notUseOcclusions) {
    ss << "noocc";
  }
  if (notUseCoplanarity) {
    ss << "_nocop";
  }
  return ss.str();
}

std::string PanoramixOptions::identityOfImage(const std::string &impath) const {
  return impath + algorithmOptionsTag();
}

void PanoramixOptions::print() const {
  std::cout << "##############################" << std::endl;
  std::cout << " useWallPrior = " << useWallPrior << std::endl;
  std::cout << " usePrincipleDirectionPrior = " << usePrincipleDirectionPrior
            << std::endl;
  std::cout << " useGeometricContextPrior = " << useGeometricContextPrior
            << std::endl;
  std::cout << " useGTOcclusions = " << useGTOcclusions << std::endl;
  std::cout << " looseLinesSecondTime = " << looseLinesSecondTime << std::endl;
  std::cout << " looseSegsSecondTime = " << looseSegsSecondTime << std::endl;
  std::cout << " restrictSegsSecondTime = " << restrictSegsSecondTime
            << std::endl;
  std::cout << " notUseOcclusions = " << notUseOcclusions << std::endl;
  std::cout << " notUseCoplanarity = " << notUseCoplanarity << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << " refresh_preparation = " << refresh_preparation << std::endl;
  std::cout << " refresh_mg_init = " << refresh_mg_init << std::endl;
  std::cout << " refresh_line2leftRightSegs = " << refresh_line2leftRightSegs
            << std::endl;
  std::cout << " refresh_mg_oriented = " << refresh_mg_oriented << std::endl;
  std::cout << " refresh_lsw = " << refresh_lsw << std::endl;
  std::cout << " refresh_mg_occdetected = " << refresh_mg_occdetected
            << std::endl;
  std::cout << " refresh_mg_reconstructed = " << refresh_mg_reconstructed
            << std::endl;
  std::cout << "##############################" << std::endl;
}

PanoramixReport::PanoramixReport() {
  time_preparation = -1;
  time_mg_init = -1;
  time_line2leftRightSegs = -1;
  time_mg_oriented = -1;
  time_lsw = -1;
  time_mg_occdetected = -1;
  time_mg_reconstructed = -1;
  succeeded = false;
}

void PanoramixReport::print() const {
  std::cout << "##############################" << std::endl;
  std::cout << " time_preparation = " << time_preparation << std::endl;
  std::cout << " time_mg_init = " << time_mg_init << std::endl;
  std::cout << " time_line2leftRightSegs = " << time_line2leftRightSegs
            << std::endl;
  std::cout << " time_mg_oriented = " << time_mg_oriented << std::endl;
  std::cout << " time_lsw = " << time_lsw << std::endl;
  std::cout << " time_mg_occdetected = " << time_mg_occdetected << std::endl;
  std::cout << " time_mg_reconstructed = " << time_mg_reconstructed
            << std::endl;
  std::cout << "##############################" << std::endl;
}

static const double thetaTiny = DegreesToRadians(2);
static const double thetaMid = DegreesToRadians(5);
static const double thetaLarge = DegreesToRadians(15);

PanoramixReport RunPanoramix(const PILayoutAnnotation &anno,
                             const PanoramixOptions &options,
                             misc::Matlab &matlab, bool showGUI,
                             bool writeToFile) {

  PanoramixReport report;
#define START_TIME_RECORD(name)                                                \
  auto start_##name = std::chrono::system_clock::now()

#define STOP_TIME_RECORD(name)                                                 \
  report.time_##name = ElapsedInMS(start_##name);                              \
  std::cout << "refresh_" #name " time cost: " << report.time_##name << "ms"   \
            << std::endl

  options.print();
  const auto identity = options.identityOfImage(anno.impath);
  misc::SaveCache(identity, "options", options);
  misc::SaveCache(identity, "report", report);

  auto image = anno.rectifiedImage.clone();
  ResizeToHeight(image, 700);

  /// prepare things!
  View<PanoramicCamera, Image3ub> view;
  std::vector<PerspectiveCamera> cams;
  std::vector<std::vector<Classified<Line2>>> rawLine2s;
  std::vector<Classified<Line3>> line3s;
  std::vector<Vec3> vps;
  int vertVPId;
  Imagei segs;
  int nsegs;

  if (options.refresh_preparation ||
      !misc::LoadCache(identity, "preparation", view, cams, rawLine2s, line3s,
                       vps, vertVPId, segs, nsegs)) {
    START_TIME_RECORD(preparation);

    view = CreatePanoramicView(image);

    // collect lines in each view
    cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows,
                                   image.rows * 0.4);
    std::vector<Line3> rawLine3s;
    rawLine2s.resize(cams.size());
    for (int i = 0; i < cams.size(); i++) {
      auto pim = view.sampled(cams[i]).image;
      LineSegmentExtractor lineExtractor;
      lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
      auto ls = lineExtractor(pim); // use pyramid
      rawLine2s[i] = ClassifyEachAs(ls, -1);
      for (auto &l : ls) {
        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                               normalize(cams[i].toSpace(l.second)));
      }
    }
    rawLine3s = MergeLines(rawLine3s, DegreesToRadians(3), DegreesToRadians(5));

    // estimate vp
    line3s = ClassifyEachAs(rawLine3s, -1);
    vps = EstimateVanishingPointsAndClassifyLines(line3s, nullptr, true);
    vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

    if (showGUI) {
      gui::ColorTable ctable = gui::RGBGreys;
      for (int i = 0; i < cams.size(); i++) {
        auto &cam = cams[i];
        std::vector<Classified<Line2>> lines;
        for (auto &l3 : line3s) {
          if (!cam.isVisibleOnScreen(l3.component.first) ||
              !cam.isVisibleOnScreen(l3.component.second)) {
            continue;
          }
          auto p1 = cam.toScreen(l3.component.first);
          auto p2 = cam.toScreen(l3.component.second);
          lines.push_back(ClassifyAs(Line2(p1, p2), l3.claz));
        }
        auto pim = view.sampled(cams[i]).image;
        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines).show();
      }
    }

    // estimate segs
    nsegs = SegmentationForPIGraph(view, line3s, segs, DegreesToRadians(1));
    RemoveThinRegionInSegmentation(segs, 1, true);
    RemoveEmbededRegionsInSegmentation(segs, true);
    nsegs = DensifySegmentation(segs, true);
    assert(IsDenseSegmentation(segs));

    if (showGUI) {
      auto ctable = gui::CreateGreyColorTableWithSize(nsegs);
      ctable.randomize();
      gui::ColorTable rgb = gui::RGBGreys;
      auto canvas = gui::MakeCanvas(view.image).alpha(0.9).add(ctable(segs));
      for (auto &l : line3s) {
        static const double sampleAngle = M_PI / 100.0;
        auto &line = l.component;
        double spanAngle = AngleBetweenDirected(line.first, line.second);
        std::vector<Point2> ps;
        ps.reserve(spanAngle / sampleAngle);
        for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
          Vec3 dir = RotateDirection(line.first, line.second, angle);
          ps.push_back(view.camera.toScreen(dir));
        }
        for (int i = 1; i < ps.size(); i++) {
          auto &p1 = ps[i - 1];
          auto &p2 = ps[i];
          if (Distance(p1, p2) >= view.image.cols / 2) {
            continue;
          }
          canvas.thickness(2);
          canvas.colorTable(rgb).add(gui::ClassifyAs(Line2(p1, p2), l.claz));
        }
      }
      canvas.show();
    }

    STOP_TIME_RECORD(preparation);

    // save
    misc::SaveCache(identity, "preparation", view, cams, rawLine2s, line3s, vps,
                    vertVPId, segs, nsegs);
  }

  //if (showGUI) {
  //  auto ctable = gui::CreateGreyColorTableWithSize(nsegs);
  //  ctable.randomize();
  //  gui::ColorTable rgb = gui::RGBGreys;
  //  rgb.exceptionalColor() = gui::Black;
  //  auto canvas = gui::MakeCanvas(view.image).alpha(0.9);
  //  for (auto &l : line3s) {
  //    static const double sampleAngle = M_PI / 100.0;
  //    auto &line = l.component;
  //    double spanAngle = AngleBetweenDirected(line.first, line.second);
  //    std::vector<Point2> ps;
  //    ps.reserve(spanAngle / sampleAngle);
  //    for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
  //      Vec3 dir = RotateDirection(line.first, line.second, angle);
  //      ps.push_back(view.camera.toScreen(dir));
  //    }
  //    for (int i = 1; i < ps.size(); i++) {
  //      auto &p1 = ps[i - 1];
  //      auto &p2 = ps[i];
  //      if (Distance(p1, p2) >= view.image.cols / 2) {
  //        continue;
  //      }
  //      canvas.thickness(2);
  //      canvas.colorTable(rgb).add(gui::ClassifyAs(Line2(p1, p2), -1));
  //    }
  //  }
  //  canvas.show();
  //}

  // gc !!!!
  std::vector<PerspectiveCamera> hcams;
  std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
  Image5d gc;
  static const int hcamNum = 16;
  static const Sizei hcamScreenSize(500, 500);
  // static const Sizei hcamScreenSize(500, 700);
  static const int hcamFocal = 200;
  std::string hcamsgcsFileName;
  {
    std::stringstream ss;
    ss << "hcamsgcs_" << hcamNum << "_" << hcamScreenSize.width << "_"
       << hcamScreenSize.height << "_" << hcamFocal;
    hcamsgcsFileName = ss.str();
  }
  if (0 || !misc::LoadCache(anno.impath, hcamsgcsFileName, hcams, gcs)) {
    // extract gcs
    hcams = CreateHorizontalPerspectiveCameras(
        view.camera, hcamNum, hcamScreenSize.width, hcamScreenSize.height,
        hcamFocal);
    gcs.resize(hcams.size());
    for (int i = 0; i < hcams.size(); i++) {
      auto pim = view.sampled(hcams[i]);
      auto pgc = ComputeIndoorGeometricContextHedau(matlab, pim.image);
      gcs[i].component.camera = hcams[i];
      gcs[i].component.image = pgc;
      gcs[i].score = abs(
          1.0 - normalize(hcams[i].forward()).dot(normalize(view.camera.up())));
    }
    misc::SaveCache(anno.impath, hcamsgcsFileName, hcams, gcs);
  }
  std::string gcmergedFileName;
  {
    std::stringstream ss;
    ss << "gc_" << hcamNum << "_" << hcamScreenSize.width << "_"
       << hcamScreenSize.height << "_" << hcamFocal;
    gcmergedFileName = ss.str();
  }
  if (0 || !misc::LoadCache(anno.impath, gcmergedFileName, gc)) {
    gc = Combine(view.camera, gcs).image;
    misc::SaveCache(anno.impath, gcmergedFileName, gc);
  }

  if (true) {
    std::vector<Imaged> gcChannels;
    cv::split(gc, gcChannels);
    auto gc3d = ConvertToImage3d(gc);
    //cv::cvtColor(gc3d, gc3d, CV_RGB2BGR);
    cv::imwrite("C:\\Users\\YANGHAO\\Pictures\\55_gc3d.png", gc3d * 255);
    gui::AsCanvas(gc3d).show(1, "gc");
  }

  // build pigraph!
  PIGraph<PanoramicCamera> mg;
  if (options.refresh_mg_init || !misc::LoadCache(identity, "mg_init", mg)) {
    std::cout << "########## refreshing mg init ###########" << std::endl;
    START_TIME_RECORD(mg_init);
    mg = BuildPIGraph(view, vps, vertVPId, segs, line3s, DegreesToRadians(1),
                      DegreesToRadians(1), DegreesToRadians(1),
                      ///!!!0.04,
                      thetaTiny, thetaLarge, thetaTiny);
    STOP_TIME_RECORD(mg_init);
    misc::SaveCache(identity, "mg_init", mg);
  }

  std::vector<std::array<std::set<int>, 2>> line2leftRightSegs;
  // static const double angleDistForSegLineNeighborhood = DegreesToRadians(5);
  if (options.refresh_line2leftRightSegs ||
      !misc::LoadCache(identity, "line2leftRightSegs", line2leftRightSegs)) {
    std::cout << "########## refreshing line2leftRightSegs ###########"
              << std::endl;
    START_TIME_RECORD(line2leftRightSegs);
    line2leftRightSegs = CollectSegsNearLines(mg, thetaMid * 2);
    STOP_TIME_RECORD(line2leftRightSegs);
    misc::SaveCache(identity, "line2leftRightSegs", line2leftRightSegs);
  }

  const auto printPIGraph = [&mg, &identity](int delay,
                                             const std::string &saveAs) {
    static const gui::ColorTable randColors =
        gui::CreateRandomColorTableWithSize(mg.nsegs);
    static gui::ColorTable rgbColors = gui::RGBGreys;
    rgbColors.exceptionalColor() = gui::Gray;
    auto pim = PrintPIGraph2(
        mg,
        [&mg](int seg, Pixel pos) -> gui::Color {
          static const gui::ColorTable ctable =
              gui::ColorTableDescriptor::RGBGreys;
          auto &c = mg.seg2control[seg];
          if (!c.used) {
            return gui::Black;
          }
          /*if (c.orientationClaz == 0) {
              return gui::Red;
          }*/
          if (c.orientationClaz != -1) {
            return ctable[c.orientationClaz].blendWith(gui::White, 0.3);
          }
          if (c.orientationNotClaz != -1) {
            static const int w = 10;
            if (IsBetween((pos.x + pos.y) % w, 0, w / 2)) {
              return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
            } else {
              return gui::White;
            }
          }
          return gui::White;
        },
        [&mg](int lp) { return gui::Transparent; },
        [&mg](int bp) -> gui::Color { return gui::Black; }, 1, 0);
    auto canvas = gui::AsCanvas(pim);
    for (auto &l : mg.lines) {
      static const double sampleAngle = M_PI / 100.0;
      auto &line = l.component;
      int claz = l.claz;
      if (claz >= mg.vps.size()) {
        claz = -1;
      }
      double spanAngle = AngleBetweenDirected(line.first, line.second);
      std::vector<Point2> ps;
      ps.reserve(spanAngle / sampleAngle);
      for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
        Vec3 dir = RotateDirection(line.first, line.second, angle);
        ps.push_back(mg.view.camera.toScreen(dir));
      }
      for (int i = 1; i < ps.size(); i++) {
        auto p1 = ToPixel(ps[i - 1]);
        auto p2 = ToPixel(ps[i]);
        if (Distance(p1, p2) >= mg.view.image.cols / 2) {
          continue;
        }
        gui::Color color = rgbColors[claz];
        cv::clipLine(cv::Rect(0, 0, canvas.image().cols, canvas.image().rows),
                     p1, p2);
        cv::line(canvas.image(), p1, p2, (cv::Scalar)color / 255.0, 2);
      }
    }
    canvas.show(delay, "pi graph");
    if (saveAs != "") {
      cv::imwrite(saveAs, Image3ub(canvas.image() * 255));
    }
    return canvas.image();
  };

  const std::string folder =
      "D:\\Panoramix\\Panorama\\images\\" + misc::Tagify(identity) + "\\";
  misc::MakeDir(folder);

  if (writeToFile) {
    cv::imwrite(folder + "im.png", anno.view.image);
  }

  if (false) {
    auto backup = mg;
    AttachPrincipleDirectionConstraints(mg);
    printPIGraph(0, folder + "principledirections.png");
    mg = backup;
    AttachWallConstraints(mg, M_PI / 60.0);
    printPIGraph(0, folder + "wall.png");
    mg = backup;
    AttachGCConstraints(mg, gc, 0.7, 0.7, true);
    printPIGraph(0, folder + "gc.png");
    mg = backup;
  }

  // attach orientation constraints
  if (options.refresh_mg_oriented ||
      !misc::LoadCache(identity, "mg_oriented", mg)) {
    std::cout << "########## refreshing mg oriented ###########" << std::endl;
    START_TIME_RECORD(mg_oriented);
    if (options.usePrincipleDirectionPrior) {
      AttachPrincipleDirectionConstraints(mg);
    }
    if (options.useWallPrior) {
      ///!!!AttachWallConstraints(mg, M_PI / 60.0);
      AttachWallConstraints(mg, thetaTiny);
    }
    if (options.useGeometricContextPrior) {
      AttachGCConstraints(mg, gc, 0.7, 0.7, true);
    }
    STOP_TIME_RECORD(mg_oriented);
    misc::SaveCache(identity, "mg_oriented", mg);
  }

  if (false) {
    Image3ub lsim = printPIGraph(0, "") * 255;
    ReverseRows(lsim);
    gui::SceneBuilder sb;
    gui::ResourceStore::set("tex", lsim);
    Sphere3 sphere;
    sphere.center = Origin();
    sphere.radius = 1.0;
    sb.begin(sphere)
        .shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama)
        .resource("tex")
        .end();
    sb.show(true, true);
  }

  // detect occlusions
  std::vector<LineSidingWeight> lsw;
  if (options.refresh_lsw || !misc::LoadCache(identity, "lsw", lsw)) {
    std::cout << "########## refreshing lsw ###########" << std::endl;
    START_TIME_RECORD(lsw);
    if (options.notUseOcclusions) {
      lsw.resize(mg.nlines(), LineSidingWeight{0.5, 0.5});
    } else {
      if (!options.useGTOcclusions) {
        lsw = ComputeLinesSidingWeights2(mg, DegreesToRadians(3), 0.2, 0.1,
                                         thetaMid);
      } else {
        lsw = ComputeLinesSidingWeightsFromAnnotation(
            mg, anno, DegreesToRadians(0.5), DegreesToRadians(8), 0.6);
      }
    }
    STOP_TIME_RECORD(lsw);
    misc::SaveCache(identity, "lsw", lsw);
  }

  if (showGUI) {
    // printPIGraph(0, folder + misc::NameOfFile(anno.impath) +
    // options.algorithmOptionsTag() + ".lines_segs.png");
  }
  const auto drawLine = [&mg](Image3f &pim, const Line3 &line,
                              const std::string &text, const gui::Color &color,
                              bool withTeeth, int linewidth, double stepAngle) {
    double angle = AngleBetweenDirected(line.first, line.second);
    std::vector<Pixel> ps;
    for (double a = 0.0; a <= angle; a += stepAngle) {
      ps.push_back(ToPixel(mg.view.camera.toScreen(
          RotateDirection(line.first, line.second, a))));
    }
    for (int i = 1; i < ps.size(); i++) {
      auto p1 = ps[i - 1];
      auto p2 = ps[i];
      if (Distance(p1, p2) >= pim.cols / 2) {
        continue;
      }
      cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
      cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, linewidth);
      if (withTeeth) {
        auto teethp =
            ToPixel(RightPerpendicularDirectiion(ecast<double>(p2 - p1))) +
            ToPixel((ecast<double>(p1) + ecast<double>(p2)) / 2.0);
        std::vector<Pixel> triangle = {p1, teethp, p2};
        cv::fillConvexPoly(pim, triangle, (cv::Scalar)color / 255.0);
      }
    }
    if (!text.empty()) {
      cv::putText(pim, text, ps.back() + Pixel(5, 0), 1, 0.7, color);
    }
  };
  if (showGUI) {
    Image3f pim1;
    {
      static gui::ColorTable ctable = gui::RGBGreys;
      ctable.exceptionalColor() = gui::Gray;
      Image3f pim = mg.view.image.clone() / 255.0;
      pim = pim * 0.3 + 0.7;
      for (int line = 0; line < mg.nlines(); line++) {
        auto &ws = lsw[line];
        auto l = mg.lines[line].component;
        drawLine(pim, l, "", ctable[mg.lines[line].claz], false, 5, 0.005);
      }
      if (writeToFile) {
        cv::imwrite(folder + "lines.png", Image3ub(pim * 255));
      }
      pim1 = pim;
    }
    Image3f pim2;
    {
      static const gui::ColorTable randColors =
          gui::CreateRandomColorTableWithSize(mg.nsegs);
      auto pim = PrintPIGraph2(
          mg,
          [&mg](int seg, Pixel pos) -> gui::Color {
            static const gui::ColorTable ctable =
                gui::ColorTableDescriptor::RGBGreys;
            auto &c = mg.seg2control[seg];
            if (!c.used) {
              return gui::White; ////
            }
            if (c.orientationClaz != -1) {
              return ctable[c.orientationClaz].blendWith(gui::White, 0.3);
            }
            if (c.orientationNotClaz != -1) {
              static const int w = 10;
              if (IsBetween((pos.x + pos.y) % w, 0, w / 2)) {
                return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
              } else {
                return gui::White;
              }
            }
            return gui::White;
          },
          [&mg](int lp) { return gui::Transparent; },
          [&mg](int bp) -> gui::Color { return gui::Black; }, 0, 0);
      // draw boundaries in a better way
      {
        Imageb bndMask(mg.segs.size(), false);
        for (auto it = mg.segs.begin(); it != mg.segs.end(); ++it) {
          auto p = it.pos();
          for (auto p2 : {Pixel(p.x + 1, p.y), Pixel(p.x, p.y + 1),
                          Pixel(p.x + 1, p.y + 1), Pixel(p.x - 1, p.y + 1)}) {
            if (!IsBetween(p2.y, 0, mg.segs.rows)) {
              continue;
            }
            p2.x = WrapBetween(p2.x, 0, mg.segs.cols);
            if (*it != mg.segs(p2)) {
              bndMask(p) = true;
            }
          }
        }
        cv::Mat element = cv::getStructuringElement(
            cv::MORPH_ELLIPSE, cv::Size(2 * 1 + 1, 2 * 1 + 1), cv::Point(1, 1));
        cv::dilate(bndMask, bndMask, element);
        for (auto it = bndMask.begin(); it != bndMask.end(); ++it) {
          if (*it) {
            pim(it.pos()) = Vec3f();
          }
        }
      }
      if (writeToFile) {
        cv::imwrite(folder + "segs.png", Image3ub(pim * 255));
      }
      Image3f im = view.image.clone() / 255.0f;
      for (int line = 0; line < mg.nlines(); line++) {
        auto &ws = lsw[line];
        auto l = mg.lines[line].component;
        if (!ws.isOcclusion()) {
          continue;
        } else if (ws.onlyConnectLeft()) {
          drawLine(pim, l, "", gui::Blue, true, 1, 0.04);
          drawLine(im, l, "", gui::Blue, true, 1, 0.04);
        } else if (ws.onlyConnectRight()) {
          drawLine(pim, l.reversed(), "", gui::Blue, true, 1, 0.04);
          drawLine(im, l.reversed(), "", gui::Blue, true, 1, 0.04);
        } else {
          drawLine(pim, l, "", gui::Gray, false, 2, 0.005);
          drawLine(im, l, "", gui::Gray, false, 2, 0.005);
        }
      }
      if (writeToFile) {
        cv::imwrite(folder + "oriented_segs_occ.png", Image3ub(pim * 255));
        cv::imwrite(folder + "occ.png", Image3ub(im * 255));
      }
      pim2 = pim;
    }

    Image3f bar(3, pim1.cols, Vec3f(1, 1, 1));
    Image3f pim;
    cv::vconcat(std::vector<Image3f>{pim1, bar, pim2}, pim);
    if (writeToFile) {
      cv::imwrite(folder + "combinedinfo.png", Image3ub(pim * 255));
    }
  }

  if (options.refresh_mg_occdetected ||
      !misc::LoadCache(identity, "mg_occdetected", mg)) {
    std::cout << "########## refreshing mg occdetected ###########"
              << std::endl;
    START_TIME_RECORD(mg_occdetected);
    ApplyLinesSidingWeights(mg, lsw, line2leftRightSegs, true);
    if (anno.extendedOnTop && !anno.topIsPlane) {
      DisableTopSeg(mg);
    }
    if (anno.extendedOnBottom && !anno.bottomIsPlane) {
      DisableBottomSeg(mg);
    }
    STOP_TIME_RECORD(mg_occdetected);
    misc::SaveCache(identity, "mg_occdetected", mg);
  }

  PIConstraintGraph cg;
  PICGDeterminablePart dp;
  if (options.refresh_mg_reconstructed ||
      !misc::LoadCache(identity, "mg_reconstructed", mg, cg, dp)) {
    std::cout << "########## refreshing mg reconstructed ###########"
              << std::endl;
    START_TIME_RECORD(mg_reconstructed);
    cg = BuildPIConstraintGraph(mg, DegreesToRadians(1), 0.01);

    // bool hasSecondTime = options.looseLinesSecondTime ||
    // options.looseSegsSecondTime || options.restrictSegsSecondTime;
    // for (int i = 0; i < (hasSecondTime ? 2 : 1); i++) {
    dp = LocateDeterminablePart(cg, DegreesToRadians(3), false);
    auto start = std::chrono::system_clock::now();
    double energy = Solve(dp, cg, matlab, 5, 1e6, !options.notUseCoplanarity);
    report.time_solve_lp = ElapsedInMS(start);
    if (IsInfOrNaN(energy)) {
      std::cout << "solve failed" << std::endl;
      return report;
    }
    //}
    STOP_TIME_RECORD(mg_reconstructed);
    misc::SaveCache(identity, "mg_reconstructed", mg, cg, dp);
  }

  // if (options.refresh_mg_reconstructed || !misc::LoadCache(identity,
  // "mg_reconstructed_connectall", mg, cg, dp)) {
  //    std::cout << "########## refreshing mg reconstructed connectall
  //    ###########" << std::endl;
  //    //START_TIME_RECORD(mg_reconstructed);
  //    cg = BuildPIConstraintGraph(mg, DegreesToRadians(1), 0.01);

  //    //bool hasSecondTime = options.looseLinesSecondTime ||
  //    options.looseSegsSecondTime || options.restrictSegsSecondTime;
  //    //for (int i = 0; i < (hasSecondTime ? 2 : 1); i++) {
  //        dp = LocateDeterminablePart(cg, DegreesToRadians(3), true);
  //        auto start = std::chrono::system_clock::now();
  //        double energy = Solve(dp, cg, matlab, 5, 1e6,
  //        !options.notUseCoplanarity);
  //        report.time_solve_lp = ElapsedInMS(start);
  //        if (IsInfOrNaN(energy)) {
  //            std::cout << "solve failed" << std::endl;
  //            return report;
  //        }
  //    //}
  //    //STOP_TIME_RECORD(mg_reconstructed);
  //    misc::SaveCache(identity, "mg_reconstructed_connectall", mg, cg, dp);
  //}

  if (false) {
    // if (options.looseLinesSecondTime) {
    if (mg.line2used.empty()) {
      mg.line2used.resize(mg.nlines(), true);
    }
    // DisorientDanglingLines3(dp, cg, mg, 0.01, 0.2);
    //}
    // if (options.looseSegsSecondTime) {
    DisorientDanglingSegs3(dp, cg, mg, 0.01, 0.3);
    //}
    // if (options.restrictSegsSecondTime) {
    // OverorientSkewSegs(dp, cg, mg, DegreesToRadians(3), DegreesToRadians(60),
    // 0.2);
    //}
    cg = BuildPIConstraintGraph(mg, DegreesToRadians(1), 0.01);
    dp = LocateDeterminablePart(cg, DegreesToRadians(3), true);

    auto start = std::chrono::system_clock::now();
    double energy = Solve(dp, cg, matlab, 5, 1e6, !options.notUseCoplanarity);
    report.time_solve_lp = ElapsedInMS(start);

    if (IsInfOrNaN(energy)) {
      std::cout << "solve failed" << std::endl;
      return report;
    }
    misc::SaveCache(identity, "mg_reconstructed2", mg, cg, dp);
  }

  if (showGUI) {
    std::pair<double, double> validRange;
    Imaged depthMap = SurfaceDepthMap(mg.view.camera, dp, cg, mg, true);
    std::vector<double> depthsArray(depthMap.begin(), depthMap.end());
    std::sort(depthsArray.begin(), depthsArray.end());
    double maxDepth = depthsArray.back();
    double minDepth = depthsArray.front();
    Imagei depthMapDisc2 = (depthMap - minDepth) / (maxDepth - minDepth) * 255;
    Imagei depthMapDisc3 = depthMap / maxDepth * 255;
    for (int &d : depthMapDisc3) {
      if (d > 255) {
        d = 255;
      }
    }
    auto jetctable = gui::CreateJetColorTableWithSize(256, gui::Black);
    auto greyctable = gui::CreateGreyColorTableWithSize(256, gui::Black);
    auto mask = GuessMask(anno);
    ResizeToHeight(mask, mg.view.image.rows);
    if (writeToFile) {
      // apply mask
      auto greyDepth = greyctable(depthMapDisc3);
      auto jetDepth = jetctable(depthMapDisc3);
      for (auto it = mask.begin(); it != mask.end(); ++it) {
        if (!*it) {
          greyDepth(it.pos()) = Vec3ub();
          jetDepth(it.pos()) = Vec3ub();
        }
      }
      cv::imwrite(folder + "depth2.png", greyDepth);
      cv::imwrite(folder + "depth.png", jetDepth);
    }

    // Imagei depthMapDisc2 =
    auto surfaceNormalMap = SurfaceNormalMap(mg.view.camera, dp, cg, mg, true);
    auto surfaceNormalMapForShow = surfaceNormalMap.clone();
    for (auto &n : surfaceNormalMapForShow) {
      auto nn = n;
      for (int i = 0; i < 3; i++) {
        n[i] = abs(nn.dot(vps[i]));
      }
      std::swap(n[0], n[2]);
    }
    if (writeToFile) {
      for (auto it = mask.begin(); it != mask.end(); ++it) {
        if (!*it) {
          surfaceNormalMapForShow(it.pos()) = Vec3();
        }
      }
      cv::imwrite(folder + "normals.png", surfaceNormalMapForShow * 255);
    }
    VisualizeReconstruction(dp, cg, mg, false,
                            [&cg, &mg](int ent) -> gui::Color {
                              auto &e = cg.entities[ent];
                              if (e.isSeg()) {
                                auto nn = normalize(
                                    e.supportingPlane.reconstructed.normal);
                                Vec3 n;
                                for (int i = 0; i < 3; i++) {
                                  n[i] = abs(nn.dot(normalize(mg.vps[i])));
                                }
                                gui::Color color = normalize(n);
                                return color;
                              } else {
                                return gui::Black;
                              }
                            },
                            nullptr, true);

    VisualizeReconstructionCompact(anno.rectifiedImage, dp, cg, mg, true);
  }

  report.succeeded = true;
  misc::SaveCache(identity, "report", report);

  return report;
}

bool GetPanoramixResult(const PILayoutAnnotation &anno,
                        const PanoramixOptions &options,
                        PIGraph<PanoramicCamera> &mg, PIConstraintGraph &cg,
                        PICGDeterminablePart &dp) {
  auto identity = options.identityOfImage(anno.impath);
  return misc::LoadCache(identity, "mg_reconstructed", mg, cg, dp);
}

std::vector<LineSidingWeight>
GetPanoramixOcclusionResult(const PILayoutAnnotation &anno,
                            const PanoramixOptions &options) {
  std::vector<LineSidingWeight> lsw;
  auto identity = options.identityOfImage(anno.impath);
  misc::LoadCache(identity, "lsw", lsw);
  return lsw;
}
