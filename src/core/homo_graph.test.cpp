#include "../core/forest.hpp"
#include "../core/homo_graph.hpp"
#include "../core/iterators.hpp"
#include "../core/utility.hpp"
#include "../gui/scene.hpp"

#include "../panoramix.unittest.hpp"

using namespace pano;
using namespace test;


TEST(ForestTest, Basic) {

  core::Forest<int> f;
  auto r = f.addRoot(0);
  auto n1 = f.add(r, 1);
  auto n2 = f.add(n1, 2);
  auto n3 = f.add(n2, 3);

  auto n12 = f.add(r, 4);

  std::vector<int> dfsResults;
  f.depthFirstSearch(r, [&dfsResults, &f](core::Forest<int>::NodeHandle nh) {
    dfsResults.push_back(f.data(nh));
    return true;
  });
  EXPECT_TRUE((dfsResults == std::vector<int>{0, 1, 2, 3, 4}));

  std::vector<int> bfsResults;
  f.breadthFirstSearch(r, [&bfsResults, &f](core::Forest<int>::NodeHandle nh) {
    bfsResults.push_back(f.data(nh));
    return true;
  });
  EXPECT_TRUE((bfsResults == std::vector<int>{0, 1, 4, 2, 3}));

  f.remove(n1);
  f.gc();

  EXPECT_EQ(f.internalNodes().size(), 2);
}

TEST(GraphicalModelTest, Basic) {

  using CGraph =
      core::HomogeneousGraph<core::Dummy,
                             core::LayerConfig<core::Dummy, core::Dynamic>>;

  CGraph cgraph;
  auto c0 = cgraph.add();
  auto c1 = cgraph.add();
  auto c2 = cgraph.add();
  auto c3 = cgraph.add();

  auto cc012 = cgraph.add<1>({c0, c1, c2});
  auto cc123 = cgraph.add<1>({c1, c2, c3});
  auto cc230 = cgraph.add<1>({c2, c3, c0});
  auto cc301 = cgraph.add<1>({c3, c0, c1});

  EXPECT_EQ(4, cgraph.internalElements<0>().size());
  EXPECT_EQ(4, cgraph.internalElements<1>().size());

  for (size_t i = 0; i < cgraph.internalElements<0>().size(); i++) {
    CGraph ncgraph = cgraph;
    ncgraph.remove(core::HandleOfTypeAtLevel<CGraph, 0>(i));
    ncgraph.gc();

    EXPECT_EQ(3, ncgraph.internalElements<0>().size());
    EXPECT_EQ(1, ncgraph.internalElements<1>().size());

    int id = 0;
    for (auto c : ncgraph.internalElements<0>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
    id = 0;
    for (auto c : ncgraph.internalElements<1>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
  }

  for (size_t i = 0; i < cgraph.internalElements<0>().size(); i++) {
    CGraph ncgraph = cgraph;
    ncgraph.remove(core::HandleOfTypeAtLevel<CGraph, 0>(i));
    ncgraph.remove(core::HandleOfTypeAtLevel<CGraph, 0>(
        (i + 1) % ncgraph.internalElements<0>().size()));
    ncgraph.gc();

    EXPECT_EQ(2, ncgraph.internalElements<0>().size());
    EXPECT_EQ(0, ncgraph.internalElements<1>().size());

    int id = 0;
    for (auto c : ncgraph.internalElements<0>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
    id = 0;
    for (auto c : ncgraph.internalElements<1>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
  }

  for (size_t i = 0; i < cgraph.internalElements<1>().size(); i++) {
    CGraph ncgraph = cgraph;
    ncgraph.remove(core::HandleOfTypeAtLevel<CGraph, 1>(i));
    ncgraph.gc();

    EXPECT_EQ(4, ncgraph.internalElements<0>().size());
    EXPECT_EQ(3, ncgraph.internalElements<1>().size());

    int id = 0;
    for (auto c : ncgraph.internalElements<0>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
    id = 0;
    for (auto c : ncgraph.internalElements<1>()) {
      EXPECT_EQ(id++, c.topo.hd.id);
    }
  }

  cgraph.remove(c0);
  cgraph.remove(c1);
  cgraph.gc();

  EXPECT_EQ(2, cgraph.internalElements<0>().size());
  EXPECT_EQ(0, cgraph.internalElements<1>().size());

  int id = 0;
  for (auto c : cgraph.internalElements<0>()) {
    EXPECT_EQ(id++, c.topo.hd.id);
  }
  id = 0;
  for (auto c : cgraph.internalElements<1>()) {
    EXPECT_EQ(id++, c.topo.hd.id);
  }
}
