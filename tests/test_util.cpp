#include "../src/core/misc.hpp"
#include "../src/core/utilities.hpp"

#include <random>
#include <list>

#include "gtest/gtest.h"

using namespace panoramix;

double randf(){
    return (std::rand() % 100000) / 100000.0;
}

TEST(UtilTest, Distance){
    using namespace core;

    auto d = Distance(std::complex<double>(1, 2), std::complex<double>(3, 4));
    ASSERT_DOUBLE_EQ(2 * sqrt(2), d);

}

TEST(UtilTest, BoundingBox) {

    using namespace core;
    Line3 l1(Point3(0.5, 0.1, 1), Point3(1, 0.4, 0.7));
    Line3 l2(Point3(0.6, 1, 0.9), Point3(0.2, -1, 0.5));

    Line3 lines[] = { l1, l2 };
    auto box = BoundingBoxOfContainer(lines);

    ASSERT_EQ(0.2, box.minCorner[0]);
    ASSERT_EQ(1, box.maxCorner[0]);

    ASSERT_EQ(-1, box.minCorner[1]);
    ASSERT_EQ(1, box.maxCorner[1]);

    ASSERT_EQ(0.5, box.minCorner[2]);
    ASSERT_EQ(1, box.maxCorner[2]);

}


TEST(UtilTest, BoundBetween) {
    for (int i = 0; i < 10000; i++){
        double x = double(rand()) / rand() + rand();
        double a = double(rand()) / rand() + rand();
        double b = a + abs(double(rand()) / rand());
        if (a == b || std::isnan(x) || std::isnan(a) || std::isnan(b) ||
            std::isinf(x) || std::isinf(a) || std::isinf(b))
            continue;
        double xx = core::BoundBetween(x, a, b);
        ASSERT_LE(a, xx);
        ASSERT_LE(xx, b);
        ASSERT_TRUE(xx == x || xx == a || xx == b);
    }
}

TEST(UtilTest, WrapBetween) {
    for (int i = 0; i < 10000; i++){
        double x = double(rand()) / rand() + rand();
        double a = double(rand()) / rand() + rand();
        double b = a + abs(double(rand())/rand());
        if (a == b || std::isnan(x) || std::isnan(a) || std::isnan(b) || 
            std::isinf(x) || std::isinf(a) || std::isinf(b))
            continue;
        double xx = core::WrapBetween(x, a, b);
        double rem = (xx - x) / (b - a) - std::round((xx - x) / (b - a));
        if (std::isnan(rem)){
            assert(0);
        }
        EXPECT_NEAR(0, rem, 1e-5);
        EXPECT_LE(a, xx);
        EXPECT_LT(xx, b);
    }

    // int test
    int x = core::WrapBetween(0, 1, 2);
    for (int i = 0; i < 10000; i++){
        int x = rand();
        int a = rand();
        int b = a + abs(rand());
        if (a == b)
            continue;

        if (core::WrapBetween(x, a, b) == b){
            std::cout << a << b << std::endl;
        }

        EXPECT_LE(a, core::WrapBetween(x, a, b));
        EXPECT_LT(core::WrapBetween(x, a, b), b);
    }
}

TEST(UtilTest, SubscriptAndIndex) {

    auto i = core::EncodeSubscriptToIndex(core::Point<int, 2>(1, 2), core::Vec<int, 2>(2, 3));
    ASSERT_EQ(5, i);

    {
        int trueIndex = 0;
        for (int a = 0; a < 10; a++){
            for (int b = 0; b < 20; b++){
                for (int c = 0; c < 15; c++){
                    for (int d = 0; d < 9; d++){
                        int index = core::EncodeSubscriptToIndex(core::Point<int, 4>(a, b, c, d),
                            core::Vec<int, 4>(10, 20, 15, 9));
                        ASSERT_EQ(trueIndex, index);

                        auto subs = core::DecodeIndexToSubscript(index, core::Vec<int, 4>(10, 20, 15, 9));
                        ASSERT_EQ(0, core::norm(core::Point<int, 4>(a, b, c, d) - subs));

                        trueIndex++;
                    }
                }
            }
        }
    }

    {
        float trueIndex = 0;
        for (int a = 0; a < 10; a++){
            for (int b = 0; b < 20; b++){
                for (int c = 0; c < 15; c++){
                    for (int d = 0; d < 9; d++){
                        float index = core::EncodeSubscriptToIndex(core::Point<float, 4>(a, b, c, d),
                            core::Vec<float, 4>(10, 20, 15, 9));
                        ASSERT_FLOAT_EQ(trueIndex, index);

                        auto subs = core::DecodeIndexToSubscript(index, core::Vec<float, 4>(10, 20, 15, 9));
                        ASSERT_FLOAT_EQ(0, core::norm(core::Point<float, 4>(a, b, c, d) - subs));

                        trueIndex++;
                    }
                }
            }
        }
    }
}


TEST(UtilTest, AngleBetweenDirections) {
    core::Vec2 v1(1, 0), v2(1, 1);
    ASSERT_DOUBLE_EQ(M_PI_4, core::AngleBetweenDirections(v1, v2));
    ASSERT_DOUBLE_EQ(M_PI_4, core::SignedAngleBetweenDirections(v1, v2, false));
    ASSERT_DOUBLE_EQ(-M_PI_4, core::SignedAngleBetweenDirections(v1, v2, true));
    ASSERT_DOUBLE_EQ(-M_PI_4, core::SignedAngleBetweenDirections(v1, v2));
    core::Vec2 v3(-1, -1);
    ASSERT_DOUBLE_EQ(-M_PI_4 * 3, core::SignedAngleBetweenDirections(v1, v3, false));
    ASSERT_DOUBLE_EQ(M_PI_4 * 3, core::SignedAngleBetweenDirections(v1, v3, true));
    ASSERT_FLOAT_EQ(M_PI, core::AngleBetweenDirections(v2, v3));
}

TEST(UtilTest, DistanceFromPointToLine) {
    core::Line3 l;
    l.first = { 1, 0, 0 };
    l.second = { -1, 0, 0 };
    for (double x = -3; x <= 3; x += 0.5) {
        core::Point3 p(x, 1, 0);
        if (x < -1){
            ASSERT_DOUBLE_EQ(core::norm(l.second - p), core::DistanceFromPointToLine(p, l).first);
        }
        else if (x > 1){
            ASSERT_DOUBLE_EQ(core::norm(l.first - p), core::DistanceFromPointToLine(p, l).first);
        }
        else{
            ASSERT_DOUBLE_EQ(1, core::DistanceFromPointToLine(p, l).first);
        }
        ASSERT_DOUBLE_EQ(1, core::norm(core::ProjectionOfPointOnLine(p, l).position, p));
    }
}

TEST(UtilTest, DistanceBetweenTwoLines) {
    core::Line3 l1 = { { 1, 0, 0 }, { -1, 0, 0 } };
    core::Line3 l2 = { { 0, 1, 1 }, { 0, -1, 1 } };
    core::Line3 l3 = { { 0, 2, 1 }, { 0, 3, 1 } };
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l2).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l2.reversed()).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1.reversed(), l2).first);
    ASSERT_DOUBLE_EQ(core::norm(l3.first - l1.center()), core::DistanceBetweenTwoLines(l1, l3).first);
    ASSERT_DOUBLE_EQ(core::norm(l3.first - l1.center()), core::DistanceBetweenTwoLines(l1.reversed(), l3).first);
    ASSERT_DOUBLE_EQ(core::norm(l3.first - l1.center()), core::DistanceBetweenTwoLines(l1, l3.reversed()).first);

    double dd = 1;
    ASSERT_DOUBLE_EQ(dd, core::DistanceBetweenTwoLines(
        core::Line3{ { 0, 0, 0 }, { 1, 0, 0 } }, 
        core::Line3{ { 0, 0, dd }, { 1, 0, dd + 0.01 } }.reversed()).first);
    ASSERT_DOUBLE_EQ(dd, core::DistanceBetweenTwoLines(
        core::Line3{ { 0, 0, 0 }, { 1, 0, 0 } },
        core::Line3{ { 0, 0, dd }, { 1, 0, dd } }.reversed()).first);

    core::Line3 l4 = { { 1, 1, 0 }, { -1, 1, 0 } };
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l4).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l4.reversed()).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1.reversed(), l4).first);

    core::Line3 l5 = { { 2, 1, 0 }, { 3, 1, 0 } };
    ASSERT_DOUBLE_EQ(core::norm(l1.first, l5.first), core::DistanceBetweenTwoLines(l1, l5).first);
    ASSERT_DOUBLE_EQ(core::norm(l1.first, l5.first), core::DistanceBetweenTwoLines(l1.reversed(), l5).first);
    ASSERT_DOUBLE_EQ(core::norm(l1.first, l5.first), core::DistanceBetweenTwoLines(l1, l5.reversed()).first);

    core::Line3 l6 = { { 2, 1, 0 }, { -2, 1, 0 } };
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l6).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1, l6.reversed()).first);
    ASSERT_DOUBLE_EQ(1, core::DistanceBetweenTwoLines(l1.reversed(), l6).first);

    core::Line3 l7 = { { 0, 0, 1 }, { 1, 0, 1 } };
    core::Line3 l8 = { { 0, 1, 0 }, { 0, 2, 0 } };
    ASSERT_DOUBLE_EQ(sqrt(2.0), core::DistanceBetweenTwoLines(l7, l8).first);

    core::Line3 l9 = { { 1, 0, 1 }, { 0, 0, 1 } };
    core::Line3 l10 = { { 0, 2, 0 }, { 0, 1, 0 } };
    ASSERT_DOUBLE_EQ(sqrt(2.0), core::DistanceBetweenTwoLines(l7, l10).first);
    ASSERT_DOUBLE_EQ(sqrt(2.0), core::DistanceBetweenTwoLines(l9, l8).first);
    ASSERT_DOUBLE_EQ(sqrt(2.0), core::DistanceBetweenTwoLines(l9, l10).first);

    ASSERT_DOUBLE_EQ(0, core::DistanceBetweenTwoLines(l7, l7).first);
    ASSERT_DOUBLE_EQ(0, core::DistanceBetweenTwoLines(l10, l10).first);

    ASSERT_DOUBLE_EQ(0, core::DistanceBetweenTwoLines(l7, l9).first);
    ASSERT_DOUBLE_EQ(0, core::DistanceBetweenTwoLines(l8, l10).first);

    core::Line3 a = { { 0.32060601460883287, 0.92543477139591090, -0.20194619899378968 }, { 0.24497237141965944, 0.95024535429166723, -0.19241180807874725 } };
    core::Line3 b = { { 0.15085395484832950, 0.90385564866472523, -0.40035990144304773 }, { 0.096251150572768340, 0.90140252138592014, -0.42214832754912596 } };
    auto pab = core::DistanceBetweenTwoLines(a, b);

    ASSERT_LE(pab.first, core::Distance(a.first, b.first));
    ASSERT_LE(pab.first, core::Distance(a.first, b.second));
    ASSERT_LE(pab.first, core::Distance(a.second, b.first));
    ASSERT_LE(pab.first, core::Distance(a.second, b.second));

    for (int i = 0; i < 1000; i++){
        core::Line3 aa = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        core::Line3 bb = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        auto p = core::DistanceBetweenTwoLines(aa, bb);

        EXPECT_LE(p.first - 0.1, core::Distance(aa.first, bb.first));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.first, bb.second));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.second, bb.first));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.second, bb.second));
    }
}


TEST(UtilTest, CreateLinearSequence) {
    size_t n = 500;
    size_t low = -150, high = 100;
    std::vector<size_t> vs(n);
    core::CreateLinearSequence(vs.begin(), vs.end(), low, high);
    for (size_t i = 0; i < n; i++){
        EXPECT_EQ(low + i * (high - low) / n, vs[i]);
    }
}

TEST(UtilTest, RTreeWrapperLargeData) {
//void run(){
    std::list<core::Line2> lines;
    std::generate_n(std::back_inserter(lines), 100000,
        [](){
        return core::Line2{ 
            core::Point2(std::rand(), std::rand()), 
            core::Point2(std::rand(), std::rand()) 
        };
    });
    core::RTreeWrapper<core::Line2> rtree(lines.begin(), lines.end());
    EXPECT_EQ(lines.size(), rtree.size());

}


TEST(UtilTest, MergeNearNaiveCorrectness) {
    std::list<double> arr1;
    arr1.resize(1000);
    std::generate(arr1.begin(), arr1.end(), std::rand);
    std::vector<double> arr2(arr1.begin(), arr1.end());

    double thres = 100;
    std::vector<std::list<double>::iterator> gBegins1;
    core::MergeNearNaive(std::begin(arr1), std::end(arr1), std::back_inserter(gBegins1), std::false_type(), thres);
    std::vector<std::vector<double>::iterator> gBegins2;
    core::MergeNearNaive(std::begin(arr2), std::end(arr2), std::back_inserter(gBegins2), std::true_type(), thres);
    ASSERT_EQ(gBegins1.size(), gBegins2.size());
    auto i = gBegins1.begin();
    auto j = gBegins2.begin();
    for (; i != gBegins1.end(); ++i, ++j){
        EXPECT_EQ(**i, **j);
    }
    for (auto i = gBegins2.begin(); i != gBegins2.end(); ++i){
        auto inext = std::next(i);
        auto begin = *i;
        auto end = inext == gBegins2.end() ? std::end(arr2) : *inext;
        auto beginVal = *begin;
        for (auto j = begin; j != end; ++j){
            EXPECT_NEAR(*j, beginVal, thres);
        }
    }
}

TEST(UtilTest, MergeNearRTreeCorrectness) {

    std::list<double> arr1;
    arr1.resize(1000);
    std::generate(arr1.begin(), arr1.end(), std::rand);
    std::vector<double> arr2(arr1.begin(), arr1.end());

    double thres = 100;
    std::vector<std::list<double>::iterator> gBegins1;
    core::MergeNearRTree(std::begin(arr1), std::end(arr1), std::back_inserter(gBegins1), std::false_type(), thres);
    std::vector<std::vector<double>::iterator> gBegins2;
    core::MergeNearNaive(std::begin(arr2), std::end(arr2), std::back_inserter(gBegins2), std::false_type(), thres);
    ASSERT_EQ(gBegins1.size(), gBegins2.size());
    auto i = gBegins1.begin();
    auto j = gBegins2.begin();
    for (; i != gBegins1.end(); ++i, ++j){
        auto iv = **i;
        auto jv = **j;
        EXPECT_DOUBLE_EQ(0, core::Distance(**i, **j));
    }

    core::RTreeWrapper<core::Line2> lines;
    lines.insert({ { 1, 2 }, { 3, 4 } });
    lines.insert({ { 2, 3 }, { 4, 5 } });

}


TEST(UtilTest, MergeNearRTreeEfficiency) {
    std::list<core::Vec4> arr1;
    std::srand(0);
    arr1.resize(5000);
    std::generate(arr1.begin(), arr1.end(), [](){ return core::Vec4(std::rand() % 100, std::rand() % 100, std::rand() % 100, std::rand() % 100); });
    double thres = 2;
    std::vector<decltype(arr1.begin())> gBegins;
    gBegins.reserve(arr1.size() / 2);
    core::MergeNearRTree(std::begin(arr1), std::end(arr1), std::back_inserter(gBegins), std::false_type(), thres);
    std::cout << gBegins.size() << std::endl;
}

TEST(UtilTest, MergeNearNaiveEfficiency) {
    std::list<core::Vec4> arr1;
    std::srand(0);
    arr1.resize(5000);
    std::generate(arr1.begin(), arr1.end(), [](){ return core::Vec4(std::rand()%100, std::rand()%100, std::rand()%100, std::rand()%100); });
    double thres = 2;
    std::vector<decltype(arr1.begin())> gBegins;
    gBegins.reserve(arr1.size() / 2);
    core::MergeNearNaive(std::begin(arr1), std::end(arr1), std::back_inserter(gBegins), std::false_type(), thres);
    std::cout << gBegins.size() << std::endl;
}






TEST(UtilTest, MinimumSpanningTree) {   
    std::vector<int> verts = { 0, 1, 2, 3, 4, 5};
    std::vector<int> edges = { 0, 1, 2, 3, 4, 5, 6, 7 };
    struct EdgeProperty { int fromv, tov; double w; };
    std::vector<EdgeProperty> edgeProperties = {
        {0, 1, 0.1},
        {1, 2, 0.2},
        {0, 2, 0.5},
        {0, 5, 0.2},
        {5, 4, 0.7},
        {2, 4, 0.3},
        {3, 4, 0.1},
        {2, 3, 0.5}
    };
    
    std::vector<int> MST;
    MST.reserve(5);
    core::MinimumSpanningTree(
        verts.begin(), verts.end(), edges.begin(), edges.end(),
        std::back_inserter(MST),
        [&edgeProperties](int e){ // get vertices of edge
            return std::make_pair(edgeProperties[e].fromv, edgeProperties[e].tov); 
        },
        [&edgeProperties](int e1, int e2){ // compare weights of edges
            return edgeProperties[e1].w < edgeProperties[e2].w; 
        }
    );
    
    std::vector<int> correctMST = { 0, 1, 3, 5, 6 };
    EXPECT_TRUE(std::is_permutation(MST.begin(), MST.end(), correctMST.begin()));
}

TEST(UtilTest, MinimumSpanningTree2) {
    std::vector<int> verts = { 0, 1, 2, 3, 4, 5, 6 };
    std::vector<int> edges = { 0, 1, 2, 3, 4, 5, 6, 7 };
    struct EdgeProperty { int fromv, tov; double w; };
    std::vector<EdgeProperty> edgeProperties = {
        { 0, 1, 0.1 },
        { 1, 2, 0.5 },
        { 2, 3, 0.2 },
        { 0, 3, 0.6 },
        { 0, 2, 0.2 },

        { 4, 5, 0.3 },
        { 4, 6, 0.8 },
        { 5, 6, 0.2 }
    };

    std::vector<int> MST;
    MST.reserve(5);
    core::MinimumSpanningTree(
        verts.begin(), verts.end(), edges.begin(), edges.end(),
        std::back_inserter(MST),
        [&edgeProperties](int e){ // get vertices of edge
        return std::make_pair(edgeProperties[e].fromv, edgeProperties[e].tov);
    },
        [&edgeProperties](int e1, int e2){ // compare weights of edges
        return edgeProperties[e1].w < edgeProperties[e2].w;
    }
    );

    std::vector<int> correctMST = { 0, 4, 2, 5, 7 };
    EXPECT_TRUE(std::is_permutation(MST.begin(), MST.end(), correctMST.begin()));
}


TEST(UtilTest, DFS) {
    std::vector<int> verts = { 0, 1, 2, 3, 4, 5, 6 };
    struct EdgeProperty { int fromv, tov; double w; };
    std::vector<EdgeProperty> edgeProperties = {
        { 0, 1, 0.1 },
        { 1, 2, 0.5 },
        { 2, 4, 0.2 },
        { 0, 4, 0.6 },
        { 0, 2, 0.2 },

        { 3, 5, 0.3 },
        { 3, 6, 0.8 },
        { 5, 6, 0.2 }
    };

    struct Data { std::shared_ptr<std::vector<int>> vntable; int vid; };
    auto compData = [](const Data & a, const Data & b) {return a.vntable == b.vntable && a.vid == b.vid; };
    auto getValue = [](const Data & cdata) -> int {return (*(cdata.vntable))[cdata.vid]; };
    auto setToNext = [](Data & cdata){cdata.vid++; };

    auto vNeighborsContainerGetter = [&verts, &edgeProperties, &compData, &getValue, &setToNext](int v) {
        std::vector<int> vns;
        for (auto & edge : edgeProperties) {
            if (edge.fromv == v)
                vns.push_back(edge.tov);
            if (edge.tov == v)
                vns.push_back(edge.fromv);
        }
        return vns;
    };

    std::vector<int> visitedVids;
    core::DepthFirstSearch(verts.begin(), verts.end(), vNeighborsContainerGetter,
        [&visitedVids](int vid){visitedVids.push_back(vid); return true; });
    std::vector<int> correctVisitedVids = {0, 1, 2, 4, 3, 5, 6};
    EXPECT_TRUE(std::equal(visitedVids.begin(), visitedVids.end(), correctVisitedVids.begin()));

    std::vector<int> ccids;
    int ccnum = core::ConnectedComponents(verts.begin(), verts.end(), vNeighborsContainerGetter,
        [&ccids](int vid, int cid){
        ccids.push_back(cid); 
    });
    std::vector<int> correctCCids = { 0, 0, 0, 0, 1, 1, 1 };

    EXPECT_EQ(2, ccnum);
    EXPECT_TRUE(std::equal(ccids.begin(), ccids.end(), correctCCids.begin()));
}



int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    //testing::FLAGS_gtest_filter = "UtilTest.DFS";
    return RUN_ALL_TESTS();
}