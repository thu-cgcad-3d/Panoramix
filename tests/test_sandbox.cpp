#include "../src/sandbox/for_expression.hpp"
#include "../src/sandbox/for_expression2.hpp"

#include "../src/deriv/data_traits_definitions.hpp"

#include "gtest/gtest.h"

template <class T>
inline T SumAll(T && t){
    return t;
}

template <class T, class ... Ts>
inline auto SumAll(T && t, Ts &&... ts)
-> decltype(t + SumAll(ts...)) {
    return t + SumAll(ts...);
}

template <class T, class ...Ts>
struct SumResult {
    using type = decltype(SumAll(std::declval<T>(), std::declval<Ts>()...));
};

template <class T, class ...Ts>
using SumResultType = typename SumResult<T, Ts...>::type;

namespace pp{
    struct K {};

    template <class A, class B> class plus{ A a; B b; };

    template <class A, class B> plus<A, B> operator + (A a, B b) {
        return plus<A, B>{a, b};
    }
}

TEST(ForExpression, TemplateTrick) {

    pp::K k1, k2, k3;
    std::cout << typeid(k1 + k2 + k3).name() << std::endl;
    std::cout << typeid(SumAll(k1, k2, k3)).name() << std::endl;
    using T = Eigen::MatrixXd;
    Eigen::MatrixXd m(2,2);

    std::cout << typeid(decltype(SumAll(m, m, m))).name() << std::endl;
    std::cout << typeid(SumResultType<T, T, T>).name() << std::endl;

    using namespace panoramix::sandbox;
    Traits<int>::print();
    Traits<double>::print();
    
    bool c = Traits<double>::c;
    Traits<Eigen::Matrix2cd>::print();
    Traits<Eigen::Array22cf>::print();
}


TEST(FunctionMatching, Cast) {

    using namespace Eigen;

    MatrixXf d;
    ArrayXXf a = ArrayXXf::Zero(2, 3);
    auto & dd = d;
    panoramix::deriv::common::Cast(MatrixXd::Ones(3, 4) * 2, dd);
    std::cout << d << std::endl;

    double aa = 5, bb;
    panoramix::deriv::common::CWiseSelect(aa, aa, bb);
    auto r = panoramix::deriv::common::CWiseSelect(d, a, aa);
    auto rr = panoramix::deriv::common::Eval(r);
    std::cout << rr << std::endl;
    auto k = panoramix::deriv::common::GeneralTranspose(a);

}

int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}