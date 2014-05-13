#include "../src/core/mesh.hpp"
#include "../src/core/mesh_maker.hpp"
#include "../src/vis/visualize3d.hpp"

#include "gtest/gtest.h"

using namespace panoramix;
using TestMesh = core::Mesh<core::Vec3>;

void foo(){
    NOT_IMPLEMENTED_YET();
}

TEST(MiscToolsTest, ConditionalIterator) {

    std::list<int> ds;
    std::generate_n(std::back_inserter(ds), 100, std::rand);
    
    auto fun = [](int dd){return dd > 50; };
    
    std::vector<int> correct;
    std::copy_if(ds.begin(), ds.end(), std::back_inserter(correct), fun);
    
    auto correctIter = correct.begin();
    for (auto d : core::MakeConditionalContainer(&ds, fun)){
        ASSERT_EQ(*correctIter, d);
        ++correctIter;
    }

}

TEST(MeshTest, Conversion) {
    using CVMesh = core::Mesh<core::Point3>;
    CVMesh mesh;
    core::MakeQuadFacedCube(mesh);
    EXPECT_EQ(8, mesh.internalVertices().size());
    EXPECT_EQ(24, mesh.internalHalfEdges().size());
    EXPECT_EQ(6, mesh.internalFaces().size());

    std::vector<core::Line3> lines;
    std::vector<core::Point3> points;
    
    mesh.remove(CVMesh::VertHandle(0));

    for (auto h : mesh.halfedges()){
        auto p1 = mesh.data(h.topo.from());
        auto p2 = mesh.data(h.topo.to());
        lines.push_back({ p1, p2 });
    }
    for (auto v : mesh.vertices()){
        points.push_back(v.data);
    }

    vis::Visualizer3D()
        << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
        << vis::manip3d::SetDefaultColor(core::ColorTag::Black)
        << lines
        << vis::manip3d::SetDefaultColor(core::ColorTag::Red)
        << vis::manip3d::SetPointSize(20)
        << points
        << vis::manip3d::SetCamera(core::PerspectiveCamera(500, 500, 500, 
            core::Vec3(-3, 0, 0), 
            core::Vec3(0.5, 0.5, 0.5)))
        << vis::manip3d::SetBackgroundColor(core::ColorTag::White)
        << vis::manip3d::AutoSetCamera
        << vis::manip3d::Show();
}

TEST(MeshTest, Tetrahedron) {

    TestMesh mesh;
    core::MakeTetrahedron(mesh);
    EXPECT_EQ(4, mesh.internalVertices().size());
    EXPECT_EQ(12, mesh.internalHalfEdges().size());
    EXPECT_EQ(4, mesh.internalFaces().size());
    
    for (size_t i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(3, nmesh.internalVertices().size());
        EXPECT_EQ(6, nmesh.internalHalfEdges().size());
        EXPECT_EQ(1, nmesh.internalFaces().size());
    }

    for (size_t i = 0; i < mesh.internalHalfEdges().size(); i++){
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::HalfHandle(i));
        nmesh.gc();

        EXPECT_EQ(4, nmesh.internalVertices().size());
        EXPECT_EQ(10, nmesh.internalHalfEdges().size());
        EXPECT_EQ(2, nmesh.internalFaces().size());
    }

    for (size_t i = 0; i < mesh.internalFaces().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::FaceHandle(i));
        nmesh.gc();

        EXPECT_EQ(4, nmesh.internalVertices().size());
        EXPECT_EQ(12, nmesh.internalHalfEdges().size());
        EXPECT_EQ(3, nmesh.internalFaces().size());
    }
    
}

TEST(MeshTest, Cube) {
    
    TestMesh mesh;
    core::MakeQuadFacedCube(mesh);
    EXPECT_EQ(8, mesh.internalVertices().size());
    EXPECT_EQ(24, mesh.internalHalfEdges().size());
    EXPECT_EQ(6, mesh.internalFaces().size());
    
    for (size_t i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(7, nmesh.internalVertices().size());
        EXPECT_EQ(18, nmesh.internalHalfEdges().size());
        EXPECT_EQ(3, nmesh.internalFaces().size());
    }

    for (size_t i = 0; i < mesh.internalHalfEdges().size(); i++){
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::HalfHandle(i));
        nmesh.gc();

        EXPECT_EQ(8, nmesh.internalVertices().size());
        EXPECT_EQ(22, nmesh.internalHalfEdges().size());
        EXPECT_EQ(4, nmesh.internalFaces().size());
    }

    for (size_t i = 0; i < mesh.internalFaces().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::FaceHandle(i));
        nmesh.gc();

        EXPECT_EQ(8, nmesh.internalVertices().size());
        EXPECT_EQ(24, nmesh.internalHalfEdges().size());
        EXPECT_EQ(5, nmesh.internalFaces().size());
    }
    
}

TEST(MeshTest, DISABLED_Sphere) {

    TestMesh mesh;
    core::MakeQuadFacedSphere(mesh, 10, 5);

}

int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}