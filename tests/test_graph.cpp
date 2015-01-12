#include "../src/core/graph.hpp"
#include "../src/core/iterators.hpp"
#include "../src/core/utilities.hpp"
#include "../src/vis/visualizers.hpp"
#include "config.hpp"

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

    vis::Visualizer viz;
    viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
    viz.installingOptions.discretizeOptions.color = vis::ColorTag::Black;
    viz.installingOptions.lineWidth = 5;
    viz.add(lines);
    viz.installingOptions.discretizeOptions.color = vis::ColorTag::Red;
    viz.installingOptions.pointSize = 20;
    viz.add(points);
    viz.renderOptions.backgroundColor = vis::ColorTag::White;
    viz.camera(core::PerspectiveCamera(500, 500, 500,
        core::Vec3(-3, 0, 0),
        core::Vec3(0.5, 0.5, 0.5)));
    viz.show();
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


TEST(ForestTest, Basic){

    core::Forest<int> f;
    auto r = f.addRoot(0);
    auto n1 = f.add(r, 1);
    auto n2 = f.add(n1, 2);
    auto n3 = f.add(n2, 3);

    auto n12 = f.add(r, 4);

    std::vector<int> dfsResults;
    f.depthFirstSearch(r, [&dfsResults, &f](core::Forest<int>::NodeHandle nh){
        dfsResults.push_back(f.data(nh)); 
        return true;
    });
    EXPECT_TRUE((dfsResults == std::vector<int>{0, 1, 2, 3, 4}));

    std::vector<int> bfsResults;
    f.breadthFirstSearch(r, [&bfsResults, &f](core::Forest<int>::NodeHandle nh){
        bfsResults.push_back(f.data(nh));
        return true;
    });
    EXPECT_TRUE((bfsResults == std::vector<int>{0, 1, 4, 2, 3}));

    f.remove(n1);
    f.gc();

    EXPECT_EQ(f.internalNodes().size(), 2);

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