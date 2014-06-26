#ifndef PANORAMIX_DERIV_DATA_TRAITS_EIGEN_HPP
#define PANORAMIX_DERIV_DATA_TRAITS_EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "glue.hpp"

#include "data_traits.hpp"
 
namespace panoramix {
    namespace deriv {

         // eigen classes
        namespace {

            template <class T>
            struct IsEigenDense : public std::is_base_of<Eigen::EigenBase<T>, T> {};

            template <class Derived>
            struct IsEigenDense<Eigen::DenseBase<Derived>> : public std::true_type {};

            template <class T>
            struct IsEigenMatrix : public std::is_base_of<Eigen::MatrixBase<T>, T> {};

            template <class Derived>
            struct IsEigenMatrix<Eigen::MatrixBase<Derived>> : public std::true_type{};

            template <class T>
            struct IsEigenArray : public std::is_base_of<Eigen::ArrayBase<T>, T> {};

            template <class Derived>
            struct IsEigenArray<Eigen::ArrayBase<Derived>> : public std::true_type{};

            template <class T>
            struct IsEigenSparseMatrix : public std::is_base_of<Eigen::SparseMatrixBase<T>, T> {};

            template <class Derived>
            struct IsEigenSparseMatrix<Eigen::SparseMatrixBase<Derived>> : public std::true_type {};

        }

        DEFINE_DATA_CATEGORY(IsEigenDense<RemoveAllType<T>>::value &&
            !IsEigenArray<RemoveAllType<T>>::value &&
            !IsEigenMatrix<RemoveAllType<T>>::value, 
            EigenDenseTag);
        DEFINE_DATA_CATEGORY(IsEigenArray<RemoveAllType<T>>::value, EigenArrayTag);
        DEFINE_DATA_CATEGORY(IsEigenMatrix<RemoveAllType<T>>::value, EigenMatrixTag);
        DEFINE_DATA_CATEGORY(IsEigenSparseMatrix<RemoveAllType<T>>::value, EigenSparseMatrixTag);

        template <class T>
        struct DataTraits<T, EigenDenseTag> {

            using StorageType = std::decay_t<decltype(std::declval<T>().eval())>;
            using ScalarType = typename StorageType::Scalar;

            static const bool shouldBeCached = std::is_same<T, StorageType>::value;

            static StorageType one() { 
                assert(false && "cannot get one from eigen dense (which is not a scalar) type!"); 
                throw ""; 
            }

            static StorageType zero() {
                assert(false && "cannot get zero from eigen dense (which is not a scalar) type!");
                throw ""; 
            }

        };

        template <class T>
        struct DataTraits<T, EigenArrayTag> : public DataTraits<T, EigenDenseTag>{
            static const RoleInProduct roleInProduct = RoleInProduct::Array;
        };

        template <class T>
        struct DataTraits<T, EigenMatrixTag> : public DataTraits<T, EigenDenseTag>{
            static const RoleInProduct roleInProduct = RoleInProduct::Matrix;
        };

        template <class T>
        struct DataTraits<T, EigenSparseMatrixTag> {
            using StorageType = std::decay_t<decltype(std::declval<T>().eval())>;
            using ScalarType = typename StorageType::Scalar;

            static const bool shouldBeCached = std::is_same<T, StorageType>::value;
            static const RoleInProduct roleInProduct = RoleInProduct::Matrix;

            static StorageType one() {
                assert(false && "cannot get one from eigen dense (which is not a scalar) type!");
                throw "";
            }

            static StorageType zero() {
                assert(false && "cannot get zero from eigen dense (which is not a scalar) type!");
                throw "";
            }
        };


        namespace common {

            // cast
            template <class D1, class D2>
            void Cast(const Eigen::MatrixBase<D1> & from, Eigen::MatrixBase<D2> & to) {
                to = from.cast<typename D2::Scalar>();
            }

            template <class D1, class D2>
            void Cast(const Eigen::ArrayBase<D1> & from, Eigen::ArrayBase<D2> & to) {
                to = from.cast<typename D2::Scalar>();
            }

            template <class D1, class D2>
            void Cast(const Eigen::MatrixBase<D1> & from, Eigen::ArrayBase<D2> & to) {
                to = from.array();
            }

            template <class D1, class D2>
            void Cast(const Eigen::ArrayBase<D1> & from, Eigen::MatrixBase<D2> & to) {
                to = from.matrix();
            }

            template <class D1, class D2>
            void Cast(const Eigen::SparseMatrixBase<D1> & from, Eigen::SparseMatrixBase<D2> & to) {
                to = from.cast<typename D2::Scalar>();
            }

            template <class D1, class D2>
            void Cast(const Eigen::SparseMatrixBase<D1> & from, Eigen::MatrixBase<D2> & to) {
                to = from.toDense();
            }



            // eval
            template <class D>
            std::decay_t<typename Eigen::DenseBase<D>::EvalReturnType> Eval(const Eigen::DenseBase<D> & d){
                return d.eval();
            }

            template <class D>
            auto Eval(const Eigen::SparseMatrixBase<D> & s) -> std::decay_t<decltype(s.eval())> {
                return s.eval();
            }
            

            // fill with scalar
            template <class D>
            std::decay_t<typename Eigen::DenseBase<D>::EvalReturnType> FillWithScalar(const Eigen::DenseBase<D> & d, typename D::Scalar s){
                auto de = d.eval();
                de.setConstant(s);
                return de;
            }

            template <class D>
            auto FillWithScalar(const Eigen::SparseMatrixBase<D> & d, typename D::Scalar s) -> std::decay_t<decltype(d.eval())> { 
                NOT_IMPLEMENTED_YET();
            }


            // general transpose
            template <class D>
            typename Eigen::DenseBase<D>::ConstTransposeReturnType GeneralTranspose(const Eigen::MatrixBase<D> & m){
                return m.transpose();
            }

            template <class D>
            D GeneralTranspose(const Eigen::ArrayBase<D> & a) {
                return a.derived();
            }

            template <class D>
            const Eigen::Transpose<D> GeneralTranspose(const Eigen::SparseMatrixBase<D> & m){
                return m.transpose();
            }


            // cwise prod
            template <class D1, class D2>
            auto CWiseProd(const Eigen::MatrixBase<D1> & a, const Eigen::MatrixBase<D2> & b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b);
            }

            template <class D1, class D2>
            auto CWiseProd(const Eigen::ArrayBase<D1> & a, const Eigen::ArrayBase<D2> & b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b);
            }

            /// TODO for sparse matrix

            // cwise select
            template <class C, class Then, class Else>
            auto CWiseSelect(const Eigen::ArrayBase<C> & cond, const Then & thenRet, const Else & elseRet)
                -> decltype((cond>0).select(thenRet, elseRet)) {
                return (cond>0).select(thenRet, elseRet);
            }

            template <class C, class Then, class Else>
            auto CWiseSelect(const Eigen::MatrixBase<C> & cond, const Then & thenRet, const Else & elseRet)
                -> decltype((cond>0).select(thenRet, elseRet)) {
                return (cond>0).select(thenRet, elseRet);
            }

            /// TODO for sparse matrix


            // sum elements
            template <class D>
            typename D::Scalar SumElements(const Eigen::DenseBase<D> & d) {
                return d.sum();
            }

            template <class D>
            typename D::Scalar SumElements(const Eigen::SparseMatrixBase<D> & d) {
                return d.sum();
            }

            // prod elements
            template <class D>
            typename D::Scalar ProdElements(const Eigen::DenseBase<D> & d) {
                return d.prod();
            }

            template <class D>
            typename D::Scalar ProdElements(const Eigen::SparseMatrixBase<D> & d) {
                return d.prod();
            }

        }


        namespace specific {

            // inverse matrix
            template <class D>
            auto InverseMatrix(const Eigen::MatrixBase<D> & m) -> decltype(m.inverse()) {
                return m.inverse();
            }   

            // dot product
            template <class D1, class D2>
            auto DotProduct(const Eigen::MatrixBase<D1> & a, const Eigen::MatrixBase<D2> & b)
                -> decltype(a.dot(b)) {
                return a.dot(b);
            }

            // cross3 prodcut
            template <class D1, class D2>
            auto Cross3Product(const Eigen::MatrixBase<D1> & a, const Eigen::MatrixBase<D2> & b)
                -> decltype(a.cross3(b)) {
                return a.cross3(b);
            }

        }


        namespace {

            template <class ScalarType>
            struct SigmoidFunctor {
                const ScalarType operator()(const ScalarType & s) const {
                    return 1.0 / (1.0 + std::exp(-s));
                }
            };

            template <class ScalarType>
            struct TanhFunctor {
                const ScalarType operator()(const ScalarType & s) const {
                    return std::tanh(s);
                }
            };

            template <class ScalarType>
            struct AcosFunctor {
                const ScalarType operator()(const ScalarType & s) const {
                    return std::acos(s);
                }
            };

        }

        namespace common {

            // sigmoid
            template <class D>
            inline auto Sigmoid(const Eigen::ArrayBase<D> & arr) -> decltype(arr.unaryExpr(panoramix::deriv::SigmoidFunctor<typename D::Scalar>())) {
                return arr.unaryExpr(panoramix::deriv::SigmoidFunctor<typename D::Scalar>());
            }

            // tanh
            template <class D>
            inline auto Tanh(const Eigen::ArrayBase<D> & arr) -> decltype(arr.unaryExpr(panoramix::deriv::TanhFunctor<typename D::Scalar>())) {
                return arr.unaryExpr(panoramix::deriv::TanhFunctor<typename D::Scalar>());
            }

            // acos
            template <class D>
            inline auto Acos(const Eigen::ArrayBase<D> & arr) -> decltype(arr.unaryExpr(panoramix::deriv::AcosFunctor<typename D::Scalar>())) {
                return arr.unaryExpr(panoramix::deriv::AcosFunctor<typename D::Scalar>());
            }

        }

    }
}

namespace std {

    template <class D>
    inline auto tanh(const Eigen::ArrayBase<D> & arr) -> decltype(panoramix::deriv::common::Tanh(arr)) {
        return panoramix::deriv::common::Tanh(arr);
    }

    template <class D>
    inline auto acos(const Eigen::ArrayBase<D> & arr) -> decltype(panoramix::deriv::common::Acos(arr)) {
        return panoramix::deriv::common::Acos(arr);
    }

}
 
#endif