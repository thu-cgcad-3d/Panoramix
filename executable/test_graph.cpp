#include "../src/core/graph.hpp"
#include "../src/core/utilities.hpp"
#include "../src/vis/visualize3d.hpp"

#include "gtest/gtest.h"
#include "test_config.hpp"

using namespace panoramix;
using namespace test;
using TestMesh = core::Mesh<core::Vec3>;


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

void VisualizeMesh(const TestMesh & mesh){
    std::vector<core::Line3> lines;
    std::vector<core::Point3> points;

    for (auto h : mesh.halfedges()){
        auto p1 = mesh.data(h.topo.from());
        auto p2 = mesh.data(h.topo.to());
        lines.push_back({ p1, p2 });
    }
    for (auto v : mesh.vertices()){
        points.push_back(v.data);
    }

    vis::Visualizer3D()
        << vis::manip3d::SetDefaultColorTable(vis::ColorTableDescriptor::RGB)
        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
        << lines
        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Red)
        << vis::manip3d::SetDefaultPointSize(20)
        << points
        << vis::manip3d::SetCamera(core::PerspectiveCamera(500, 500, 500,
        core::Vec3(-3, 0, 0),
        core::Vec3(0.5, 0.5, 0.5)))
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
        << vis::manip3d::Show();
}


TEST(MeshTest, Visualize) {
    using namespace std::placeholders;
    TestMesh mesh;
    auto vertexAdder = [&mesh](double x, double y, double z){
        return mesh.addVertex(core::Point3(x, y, z));
    };
    auto triFaceAdder = [&mesh](const TestMesh::VertHandle a, const TestMesh::VertHandle b,
        const TestMesh::VertHandle c){
        return mesh.addFace(a, b, c);
    };
    auto quadFaceAdder = [&mesh](const TestMesh::VertHandle a, const TestMesh::VertHandle b,
        const TestMesh::VertHandle c, const TestMesh::VertHandle d){
        return mesh.addFace(a, b, c, d);
    };

    core::MakeQuadFacedCube(vertexAdder, quadFaceAdder);
    EXPECT_EQ(8, mesh.internalVertices().size());
    EXPECT_EQ(24, mesh.internalHalfEdges().size());
    EXPECT_EQ(6, mesh.internalFaces().size());

    mesh.remove(TestMesh::VertHandle(0));
    mesh.gc();
    VisualizeMesh(mesh);

    mesh.clear();
    core::MakeQuadFacedSphere(vertexAdder, quadFaceAdder, 10, 20);
    VisualizeMesh(mesh);

    mesh.clear();
    core::MakeIcosahedron(vertexAdder, triFaceAdder);
    VisualizeMesh(mesh);

    mesh.clear();
    core::MakeTriMeshFromSMFFile(vertexAdder, triFaceAdder, ProjectDataDirStrings::MeshSMF + "/bones-clean.smf");
    VisualizeMesh(mesh);

    
}

TEST(MeshTest, Tetrahedron) {

    TestMesh mesh;
    auto vertexAdder = [&mesh](double x, double y, double z){
        return mesh.addVertex(core::Point3(x, y, z));
    };
    auto triFaceAdder = [&mesh](const TestMesh::VertHandle a,
        const TestMesh::VertHandle b, const TestMesh::VertHandle c){
        return mesh.addFace({ a, b, c });
    };
    core::MakeTetrahedron(vertexAdder, triFaceAdder);
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
    auto vertexAdder = [&mesh](double x, double y, double z){
        return mesh.addVertex(core::Point3(x, y, z));
    };
    auto quadFaceAdder = [&mesh](const TestMesh::VertHandle a, const TestMesh::VertHandle b,
        const TestMesh::VertHandle c, const TestMesh::VertHandle d){
        return mesh.addFace({ a, b, c, d });
    };
    core::MakeQuadFacedCube(vertexAdder, quadFaceAdder);
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


TEST(ConstraintGraphTest, Basic) {
    
    using CGraph = core::ConstraintGraph<core::Dummy, core::Dummy>;

    CGraph cgraph;
    auto c0 = cgraph.addComponent();
    auto c1 = cgraph.addComponent();
    auto c2 = cgraph.addComponent();
    auto c3 = cgraph.addComponent();

    auto cc012 = cgraph.addConstraint({ c0, c1, c2 });
    auto cc123 = cgraph.addConstraint({ c1, c2, c3 });
    auto cc230 = cgraph.addConstraint({ c2, c3, c0 });
    auto cc301 = cgraph.addConstraint({ c3, c0, c1 });

    EXPECT_EQ(4, cgraph.internalComponents().size());
    EXPECT_EQ(4, cgraph.internalConstraints().size());
    
    for (size_t i = 0; i < cgraph.internalComponents().size(); i++){
        CGraph ncgraph = cgraph;
        ncgraph.remove(CGraph::ComponentHandle(i));
        ncgraph.gc();

        EXPECT_EQ(3, ncgraph.internalComponents().size());
        EXPECT_EQ(1, ncgraph.internalConstraints().size());
    }

    for (size_t i = 0; i < cgraph.internalConstraints().size(); i++){
        CGraph ncgraph = cgraph;
        ncgraph.remove(CGraph::ConstraintHandle(i));
        ncgraph.gc();

        EXPECT_EQ(4, ncgraph.internalComponents().size());
        EXPECT_EQ(3, ncgraph.internalConstraints().size());
    }

    cgraph.remove(c0);
    cgraph.remove(c1);
    cgraph.gc();

    EXPECT_EQ(2, cgraph.internalComponents().size());
    EXPECT_EQ(0, cgraph.internalConstraints().size());

}




TEST(HeterogeneousBinaryGraph, HeterogeneousVector){

    core::HeterogeneousVector<double, int, bool> v;
    v.pushBack(1);
    v.pushBack(2.5);
    v.pushBack(true);

    ASSERT_EQ(3, v.size());
    ASSERT_EQ(1, v.size<int>());
}



TEST(GraphicalModelTest, Basic) {

    using CGraph = core::HomogeneousGraph<core::Dummy, core::LayerConfig<core::Dummy, core::Dynamic>>;

    CGraph cgraph;
    auto c0 = cgraph.add();
    auto c1 = cgraph.add();
    auto c2 = cgraph.add();
    auto c3 = cgraph.add();

    auto cc012 = cgraph.add<1>({ c0, c1, c2 });
    auto cc123 = cgraph.add<1>({ c1, c2, c3 });
    auto cc230 = cgraph.add<1>({ c2, c3, c0 });
    auto cc301 = cgraph.add<1>({ c3, c0, c1 });

    EXPECT_EQ(4, cgraph.internalElements<0>().size());
    EXPECT_EQ(4, cgraph.internalElements<1>().size());

    for (size_t i = 0; i < cgraph.internalElements<0>().size(); i++){
        CGraph ncgraph = cgraph;
        ncgraph.remove(core::HandleAtLevel<0>(i));
        ncgraph.gc();

        EXPECT_EQ(3, ncgraph.internalElements<0>().size());
        EXPECT_EQ(1, ncgraph.internalElements<1>().size());

        int id = 0;
        for (auto c : ncgraph.internalElements<0>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
        id = 0;
        for (auto c : ncgraph.internalElements<1>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
    }

    for (size_t i = 0; i < cgraph.internalElements<0>().size(); i++){
        CGraph ncgraph = cgraph;
        ncgraph.remove(core::HandleAtLevel<0>(i));
        ncgraph.remove(core::HandleAtLevel<0>((i+1) % ncgraph.internalElements<0>().size()));
        ncgraph.gc();

        EXPECT_EQ(2, ncgraph.internalElements<0>().size());
        EXPECT_EQ(0, ncgraph.internalElements<1>().size());

        int id = 0;
        for (auto c : ncgraph.internalElements<0>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
        id = 0;
        for (auto c : ncgraph.internalElements<1>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
    }

    for (size_t i = 0; i < cgraph.internalElements<1>().size(); i++){
        CGraph ncgraph = cgraph;
        ncgraph.remove(core::HandleAtLevel<1>(i));
        ncgraph.gc();

        EXPECT_EQ(4, ncgraph.internalElements<0>().size());
        EXPECT_EQ(3, ncgraph.internalElements<1>().size());

        int id = 0;
        for (auto c : ncgraph.internalElements<0>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
        id = 0;
        for (auto c : ncgraph.internalElements<1>()){
            EXPECT_EQ(id++, c.topo.hd.id);
        }
    }

    cgraph.remove(c0);
    cgraph.remove(c1);
    cgraph.gc();

    EXPECT_EQ(2, cgraph.internalElements<0>().size());
    EXPECT_EQ(0, cgraph.internalElements<1>().size());

    int id = 0;
    for (auto c : cgraph.internalElements<0>()){
        EXPECT_EQ(id++, c.topo.hd.id);
    }
    id = 0;
    for (auto c : cgraph.internalElements<1>()){
        EXPECT_EQ(id++, c.topo.hd.id);
    }
}



int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}