#include "../src/core/any.hpp"
#include "../src/core/iterators.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/core/containers.hpp"

#include <random>
#include <list>

#include "config.hpp"

using namespace panoramix;

inline double randf(){
    return (std::rand() % 100000) / 100000.0;
}


TEST(MiscTest, Any) {

    using namespace core;

    Any aInt = 1;
    Any aFloat = 2.0f;
    Any aVector = std::vector<std::array<int, 3>>(5);

    ASSERT_TRUE(aInt.is<int>());
    ASSERT_TRUE(aFloat.is<float>());
    ASSERT_TRUE((aVector.is<std::vector<std::array<int, 3>>>()));
    ASSERT_FALSE((aVector.is<std::vector<std::array<int, 2>>>()));

    int vInt = aInt;
    ASSERT_EQ(vInt, 1);
    float vFloat = aFloat;
    ASSERT_EQ(vFloat, 2.0f);
    std::vector<std::array<int, 3>> vVector = aVector;
    ASSERT_TRUE((vVector == std::vector<std::array<int, 3>>(5)));

    std::swap(aInt, aFloat);
    ASSERT_EQ((float)aInt, 2.0f);
    ASSERT_EQ((int)aFloat, 1);

    ASSERT_FALSE(aVector.null());
    aVector = nullptr;
    ASSERT_TRUE(aVector.null());

    Any something;
    ASSERT_TRUE(something.null());
    something = std::list<int>{ 1, 2, 3, 4 };
    ASSERT_FALSE(something.null());
    ASSERT_TRUE(something.is<std::list<int>>());
    something = std::set<double>{1.0, 2.0, 3.0};
    ASSERT_FALSE(something.is<std::list<int>>());
    ASSERT_TRUE(something.is<std::set<double>>());

    // see http://stackoverflow.com/questions/20165166/double-delete-in-initializer-list-vs-2013
    Any mess = std::list<Any>{
        std::make_tuple(1, 3.0, "123"),
            std::vector<Any> {
            std::array<int, 3>{{ 1, 2, 3 }},
                std::pair<double, double>{ 1.5, 7.4 }
        }
    };

    ASSERT_TRUE((mess.cast<std::list<Any>>().back()
        .cast<std::vector<Any>>().front().cast<std::array<int, 3>>()
        == std::array<int, 3>{{ 1, 2, 3 }}));
    ASSERT_TRUE((mess.cast<std::list<Any>>().back()
        .cast<std::vector<Any>>().back().cast<std::pair<double, double>>()
        == std::pair<double, double>{ 1.5, 7.4 }));


}




TEST(MiscTest, Failable){

    core::Failable<std::vector<int>> opt;
    ASSERT_TRUE(opt.null());
    opt = core::AsResult(std::vector<int>{1, 2, 3, 4});
    ASSERT_TRUE(!opt.null());

    auto data = opt.unwrap();
    ASSERT_TRUE((data == std::vector<int>{1, 2, 3}));
    ASSERT_TRUE(opt.null());

    auto opt2 = core::AsResult(std::move(data));
    ASSERT_TRUE(!opt2.null());

    opt = std::move(opt2);
    ASSERT_TRUE(!opt.null());

}


TEST(MistTest, YieldIterator) {

    std::vector<int> data(10);
    std::iota(data.begin(), data.end(), 0);
    std::copy(data.begin(), data.end(), core::MakeYield<int>([](int d){
        std::cout << d << std::endl;
    }));

}


TEST(UtilTest, HasValue){

    std::vector<core::Ratio<core::Line2, double>> hlines = {
        {
            { {1.0, 2.0}, {4.0, 5.0} },
            0.0
        },
        {
            { { 1.0, 2.0 }, { NAN, 5.0 } },
            0.0
        },
        {
            { { 1.0, 2.0 }, { 0.0, 5.0 } },
            0.0
        },
        {
            { { 1.0, 2.0 }, { 4.0, 5.0 } },
            0.0
        },
        {
            { { 1.0, 2.0 }, { 4.0, 5.0 } },
            0.0
        }
    };

    ASSERT_TRUE(core::HasValue(hlines, std::isnan<double>));
    ASSERT_FALSE(core::HasValue(hlines, std::isinf<double>));
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

    Enabled<Line3> lines2[] = { EnableAs(l1, true), EnableAs(l2, false) };
    ASSERT_TRUE(BoundingBoxOfContainer(lines2) == BoundingBox(l1));

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

        int xx = core::WrapBetween(x, a, b);
        int rem = (xx - x) % (b - a);

        EXPECT_EQ(0, rem);
        EXPECT_LE(a, xx);
        EXPECT_LT(xx, b);        
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
    ASSERT_DOUBLE_EQ(0.0, core::AngleBetweenDirections(core::Vec3(1.0, 1.9, 0.1), core::Vec3(1.0, 1.9, 0.1000000001)));
}

TEST(UtilTest, DistanceFromPointToLine) {
    core::Line3 l;
    l.first = { 1, 0, 0 };
    l.second = { -1, 0, 0 };
    auto infLine = l.infiniteLine();
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
        ASSERT_DOUBLE_EQ(1, core::norm(core::ProjectionOfPointOnLine(p, l).position - p));
        ASSERT_DOUBLE_EQ(1, core::norm(core::DistanceFromPointToLine(p, l.infiniteLine()).second - p));
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

    {
        core::Line3 a = { { 0.32060601460883287, 0.92543477139591090, -0.20194619899378968 }, { 0.24497237141965944, 0.95024535429166723, -0.19241180807874725 } };
        core::Line3 b = { { 0.15085395484832950, 0.90385564866472523, -0.40035990144304773 }, { 0.096251150572768340, 0.90140252138592014, -0.42214832754912596 } };
        auto pab = core::DistanceBetweenTwoLines(a, b);

        ASSERT_LE(pab.first, core::Distance(a.first, b.first));
        ASSERT_LE(pab.first, core::Distance(a.first, b.second));
        ASSERT_LE(pab.first, core::Distance(a.second, b.first));
        ASSERT_LE(pab.first, core::Distance(a.second, b.second));
    }
    {
        core::Line3 a = { { 0.98184, -0.120335, 0.146665 }, { 0.65886, 0.72241, -0.209827 } };
        core::Line3 b = { { 0.493696, 0.844419, 0.207651 }, { 0.245523, 0.952812, 0.178513 } };
        auto pab = core::DistanceBetweenTwoLines(a, b);

        ASSERT_LE(pab.first, core::Distance(a.first, b.first));
        ASSERT_LE(pab.first, core::Distance(a.first, b.second));
        ASSERT_LE(pab.first, core::Distance(a.second, b.first));
        ASSERT_LE(pab.first, core::Distance(a.second, b.second));
    }

    for (int i = 0; i < 1000; i++){
        core::Line3 aa = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        core::Line3 bb = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        auto p = core::DistanceBetweenTwoLines(aa, bb);
        core::Ray3 iaa = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        core::Ray3 ibb = { { randf(), randf(), randf() }, { randf(), randf(), randf() } };
        auto ip = core::DistanceBetweenTwoLines(iaa, ibb);

        EXPECT_LE(p.first - 0.1, core::Distance(aa.first, bb.first));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.first, bb.second));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.second, bb.first));
        EXPECT_LE(p.first - 0.1, core::Distance(aa.second, bb.second));

        EXPECT_LE(ip.first - 0.1, core::Distance(iaa.anchor, ibb.anchor));
    }
}


TEST(UtilTest, BarycentricCoordinatesOfLineAndPlaneUnitIntersection) {

    core::Point3 pts[] = {
        core::Point3(1, 0, 0),
        core::Point3(0, 1, 0),
        core::Point3(0, 0, 1)
    };

    auto bc1 = core::BarycentricCoordinatesOfLineAndPlaneUnitIntersection(
        core::Ray3(core::Point3(0, 0, 0), core::Point3(1, 1, 1)), pts);

    EXPECT_FLOAT_EQ(bc1[0], 1.0/3.0);
    EXPECT_FLOAT_EQ(bc1[1], 1.0/3.0);
    EXPECT_FLOAT_EQ(bc1[2], 1.0/3.0);

    for (int i = 0; i < 3; i++){
        auto bc2 = core::BarycentricCoordinatesOfLineAndPlaneUnitIntersection(
            core::Ray3(core::Point3(0, 0, 0), pts[i]), pts);
        for (int k = 0; k < 3; k++){
            EXPECT_FLOAT_EQ(bc2[k], i == k ? 1.0 : 0.0);
        }
    }

    for (int i = 0; i < 3; i++){
        core::Vec3 dir(1, 1, 1);
        dir(i) = -.5;
        auto bc3 = core::BarycentricCoordinatesOfLineAndPlaneUnitIntersection(
            core::Ray3(core::Point3(0, 0, 0), dir), pts);
        for (int k = 0; k < 3; k++){
            EXPECT_FLOAT_EQ(bc3[k], k == i ? -1.0 / 3.0 : 2.0 / 3.0);
        }
    }

}


TEST(UtilTest, EigenVectorsAndValues) {
    {
        std::vector<core::Point3> pts;
        std::generate_n(std::back_inserter(pts), 10000, [](){
            auto a = randf();
            auto b = randf();
            auto c = randf();
            return core::Point3(a, b / 1000.0, c / 1000.0);
        });
        std::generate_n(std::back_inserter(pts), 1, [](){
            auto a = randf();
            auto b = randf();
            auto c = randf();
            return core::Point3(a, b / 1000.0 + 1.0, c / 1000.0);
        });
        auto result = core::EigenVectorAndValuesFromPoints(pts);
        for (auto & r : result){
            std::cout << "eigen vector: " << r.component << "  value: " << r.score << std::endl;
        }

        std::sort(result.begin(), result.end());
        std::cout << "planarity = " << (result[1].score / result[2].score * result[1].score / result[0].score) << std::endl;
    }
    
    core::Plane3 plane(core::Point3(1, 2, 3), core::Vec3(-2, -5, -7));
    core::Vec3 x, y;
    std::tie(x, y) = core::ProposeXYDirectionsFromZDirection(plane.normal);
    std::vector<core::Point3> pts;
    std::generate_n(std::back_inserter(pts), 10000, [&](){
        auto a = randf();
        auto b = randf();
        auto c = randf();
        return a * x + b * y + c/1000.0 * plane.normal + plane.anchor;
    });
    auto result = core::EigenVectorAndValuesFromPoints(pts);
    for (auto & r : result){
        std::cout << "eigen vector: " << core::normalize(r.component) << "  value: " << r.score << std::endl;
    }
    std::sort(result.begin(), result.end());
    std::cout << "planarity = " << (result[1].score / result[2].score * result[1].score / result[0].score) << std::endl;
}

TEST(ContainerTest, RTreeWrapperLargeData) {

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


TEST(ContainerTest, VecSetAndVecMap){

    core::VecMultiSet<double, 3> s(0.1);
    int N = 5000;
    for (int i = 0; i < N; i++){
        s.insert(core::Vec3(randf(), randf(), randf()));
    }
    
    EXPECT_EQ(s.fullSize(), N);
    for (auto & vs : s){
        for (auto & v : vs){
            ASSERT_LE(core::Distance(v, vs.front()), s.influenceRange());
        }
    }

    core::VecMap<double, 3, std::vector<int>> m(0.1);
    std::vector<core::Vec3> vecs(N);
    std::generate(vecs.begin(), vecs.end(), [](){return core::Vec3(randf(), randf(), randf()); });
    for (int i = 0; i < N; i++){
        m[vecs[i]].push_back(i);
    }

    for (auto & g : m){
        for (int i : g.second){
            ASSERT_LE(core::Distance(vecs[i], vecs[g.second.front()]), m.influenceRange());
        }
    }

}


TEST(ContainerTest, MaxHeap){
    std::vector<double> data(50000);
    std::generate(data.begin(), data.end(), randf);
    std::vector<int> ids(data.size());
    std::iota(ids.begin(), ids.end(), 0);

    std::vector<core::Scored<int>> qd(data.size());
    for (int i = 0; i < data.size(); i++)
        qd[i] = core::ScoreAs(i, data[i]);
    std::priority_queue<core::Scored<int>> Q(qd.begin(), qd.end());
    core::MaxHeap<int> H(ids.begin(), ids.end(), [&data](int id){return data[id]; });

    ASSERT_EQ(Q.size(), H.size());

    int count = 0;
    while (!Q.empty()){
        ASSERT_EQ(Q.size(), H.size());
        ASSERT_EQ(Q.top().score, H.topScore());
        Q.pop();
        H.pop();

        if (count % 2 == 0){
            double v = randf();
            Q.push(core::ScoreAs(count, v));
            H.push(count, v);
        }

        count++;
    }

    core::MaxHeap<int> HH;
    int N = 5000;
    for (int i = 0; i < N; i++){
        HH.push(i, randf());
    }
    for (int i = 0; i < N * 3; i++){
        int key = i % N;
        HH.setScore(key, randf());
        int topKey = HH.top();
        // assert topKey has the highest score
        for (int j = 0; j < N; j++){
            ASSERT_LE(HH.at(j), HH.at(topKey));
        }
    }
}



TEST(AlgorithmsTest, ForeachCompatible){
    std::list<double> data = { 1.0, 2.0, 3.0, 5.0, 6.0, 7.0 };
    std::list<double> selected;
    core::ForeachCompatibleWithLastElement(data.begin(), data.end(), std::back_inserter(selected), 
        [](double a, double b){
        return abs(a - b) >= 1.5; 
    });
    std::list<double> groundTruth = { 1.0, 3.0, 5.0, 7.0 };
    ASSERT_TRUE(selected == groundTruth);
}


TEST(AlgorithmsTest, MergeNearNaiveCorrectness) {
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

TEST(AlgorithmsTest, MergeNearRTreeCorrectness) {

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


TEST(AlgorithmsTest, MergeNearRTreeEfficiency) {
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

TEST(AlgorithmsTest, MergeNearNaiveEfficiency) {
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






TEST(AlgorithmsTest, MinimumSpanningTree) {
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

TEST(AlgorithmsTest, MinimumSpanningTree2) {
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


TEST(AlgorithmsTest, DFS_CC) {
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


TEST(AlgorithmsTest, TopologicalSort){

    {
        std::vector<int> verts = { 0, 1, 2, 3, 4, 5, 6 };
        std::random_shuffle(verts.begin(), verts.end());
        struct Edge { int from, to; };
        std::vector<Edge> edges = {
                { 0, 1 },
                { 0, 2 },
                { 1, 3 },
                { 2, 4 },
                { 1, 6 },
                { 4, 5 }
        };
        std::vector<int> sortedVerts;
        core::TopologicalSort(verts.begin(), verts.end(), std::back_inserter(sortedVerts), [&edges](int vert){
            std::vector<int> predecessors;
            for (auto & e : edges){
                if (e.to == vert)
                    predecessors.push_back(e.from);
            }
            return predecessors;
        });
        for (auto & e : edges){
            auto fromPos = std::find(sortedVerts.begin(), sortedVerts.end(), e.from) - sortedVerts.begin();
            auto toPos = std::find(sortedVerts.begin(), sortedVerts.end(), e.to) - sortedVerts.begin();
            ASSERT_TRUE(fromPos <= toPos);
        }
    }

    {
        std::vector<int> verts(1000);
        std::iota(verts.begin(), verts.end(), 0);
        std::random_shuffle(verts.begin(), verts.end());
        struct Edge { int from, to; };
        std::vector<Edge> edges(1000);
        std::generate(edges.begin(), edges.end(), [&verts](){
            int v1 = rand() % verts.size();
            int v2 = rand() % verts.size();
            return v1 < v2 ? Edge{ v1, v2 } : Edge{ v2, v1 };
        });
        std::vector<int> sortedVerts;
        core::TopologicalSort(verts.begin(), verts.end(), std::back_inserter(sortedVerts), [&edges](int vert){
            std::vector<int> predecessors;
            for (auto & e : edges){
                if (e.to == vert)
                    predecessors.push_back(e.from);
            }
            return predecessors;
        });
        for (auto & e : edges){
            auto fromPos = std::find(sortedVerts.begin(), sortedVerts.end(), e.from) - sortedVerts.begin();
            auto toPos = std::find(sortedVerts.begin(), sortedVerts.end(), e.to) - sortedVerts.begin();

            ASSERT_TRUE(fromPos <= toPos);
        }
    }

}

//
//int main(int argc, char * argv[], char * envp[])
//{
//    testing::InitGoogleTest(&argc, argv);
//    testing::GTEST_FLAG(catch_exceptions) = false;
//    testing::GTEST_FLAG(throw_on_failure) = true;
//    testing::GTEST_FLAG(filter) = "UtilTest.BarycentricCoordinatesOfLineAndPlaneUnitIntersection";
//    return RUN_ALL_TESTS();
//}