#include "gtest/gtest.h"

#include "../src/gui/singleton.hpp"
#include "config.hpp"

namespace pano {
namespace test {

std::string ProjectDataDirStrings::Base = PROJECT_TEST_DATA_DIR_STR;
std::string ProjectDataDirStrings::Normal = PROJECT_TEST_DATA_DIR_STR "/normal";
std::string ProjectDataDirStrings::PanoramaIndoor =
    PROJECT_TEST_DATA_DIR_STR "/panorama/indoor";
std::string ProjectDataDirStrings::PanoramaOutdoor =
    PROJECT_TEST_DATA_DIR_STR "/panorama/outdoor";
std::string ProjectDataDirStrings::Serialization =
    PROJECT_TEST_DATA_DIR_STR "/serialization";
std::string ProjectDataDirStrings::LocalManhattan =
    PROJECT_TEST_DATA_DIR_STR "/localmanh";
std::string ProjectDataDirStrings::MeshSMF =
    PROJECT_TEST_DATA_DIR_STR "/mesh/smf";
std::string ProjectDataDirStrings::Scripts =
    PROJECT_TEST_DATA_DIR_STR "/scripts";
std::string ProjectDataDirStrings::BPTests =
    PROJECT_TEST_DATA_DIR_STR "/bptests";
}
}

int main(int argc, char *argv[], char *envp[]) {
  pano::gui::Singleton::SetCmdArgs(argc, argv, envp);
  testing::InitGoogleTest(&argc, argv);

  // testing::GTEST_FLAG(catch_exceptions) = false;
  // testing::GTEST_FLAG(throw_on_failure) = true;
  // testing::GTEST_FLAG(filter) = "MatlabEngine.*";
  // testing::GTEST_FLAG(filter) = "BasicType.Line";
  // testing::GTEST_FLAG(filter) = "BasicType.VecCast";
  // testing::GTEST_FLAG(filter) = "MiscTest.Failable";
  // testing::GTEST_FLAG(filter) = "Feature.SegmentationExtractor";
  // testing::GTEST_FLAG(filter) = "Feature.VanishingPointsDetector";
  // testing::GTEST_FLAG(filter) = "Feature.IndoorGeometricContextMatlab";
  // testing::GTEST_FLAG(filter) = "Engine.RebuildOnOnePanorama";
  // testing::GTEST_FLAG(filter) = "Matlab.SparseMatrix";
  // testing::GTEST_FLAG(filter) = "View.OrientationContext";
  // testing::GTEST_FLAG(filter) = "ConstraintGraph.Basic";
  // testing::GTEST_FLAG(filter) = "View.OrientationContextPerspective";
  // testing::GTEST_FLAG(filter) =
  // "View.SampleViews:View.LinesGraph:View.OrientationContext";
  // testing::GTEST_FLAG(filter) = "Serialization.SUNAnnotation";
  // testing::GTEST_FLAG(filter) = "BasicTest.Move";
  // testing::GTEST_FLAG(filter) = "FactorGraph.Denoise";
  // testing::GTEST_FLAG(filter) = "Feature.SegmentationBoundaryJunction";
  // testing::GTEST_FLAG(filter) = "Serialization.Mex";
  // testing::GTEST_FLAG(filter) = "Engine.Batch";
  // testing::GTEST_FLAG(filter) = "Camera.*";
  // testing::GTEST_FLAG(filter) = "Camera.UniformSphericalCamera";
  // testing::GTEST_FLAG(filter) = "UtilTest.PlaneIntersection";
   testing::GTEST_FLAG(filter) = "Scene.LayeredShape3";
  // testing::GTEST_FLAG(filter) = "Scene.Interaction";
  // testing::GTEST_FLAG(filter) = "PPattern.*";
  // testing::GTEST_FLAG(filter) = "BasicType.*";
  // testing::GTEST_FLAG(filter) = "Feature.OcclusionDetection";
  // testing::GTEST_FLAG(filter) = "ContainerTest.Dictionary";
  // testing::GTEST_FLAG(filter) = "ContainerTest.MaxHeap";
  // testing::GTEST_FLAG(filter) = "SingleView.ComputeSpatialRegionProperties";
  // testing::GTEST_FLAG(filter) = "Feature.ExtractSegmentationTopology";
  // testing::GTEST_FLAG(filter) = "Feature.RemoveThinRegionInSegmentation";
  // testing::GTEST_FLAG(filter) = "Feature.RemoveSmallRegionInSegmentation";
  // testing::GTEST_FLAG(filter) = "Feature.SegmentationExtractorInPanorama";

  return RUN_ALL_TESTS();
}
