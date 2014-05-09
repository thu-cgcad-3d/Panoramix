#ifndef PANORAMIX_DERIV_DATA_TRAITS_EIGEN_HPP
#define PANORAMIX_DERIV_DATA_TRAITS_EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

#include "glue.hpp"

#include "data_traits.hpp"
 
namespace panoramix {
    namespace deriv {

         // eigen classes
        namespace {

            template <class T>
            struct IsEigenDense
                : public std::is_base_of<Eigen::EigenBase<T>, T>
            {};

            template <class Derived>
            struct IsEigenDense<Eigen::DenseBase<Derived>> : public std::true_type{};

            template <class T>
            struct IsEigenMatrix
                : public std::is_base_of<Eigen::MatrixBase<T>, T>
            {};

            template <class Derived>
            struct IsEigenMatrix<Eigen::MatrixBase<Derived>> : public std::true_type{};

            template <class T>
            struct IsEigenArray
                : public std::is_base_of<Eigen::ArrayBase<T>, T>
            {};

            template <class Derived>
            struct IsEigenArray<Eigen::ArrayBase<Derived>> : public std::true_type{};

        }

        DEFINE_DATA_CATEGORY(IsEigenDense<RemoveAllType<T>>::value &&
            !IsEigenArray<RemoveAllType<T>>::value &&
            !IsEigenMatrix<RemoveAllType<T>>::value,
            EigenDenseTag);
        DEFINE_DATA_CATEGORY(IsEigenArray<RemoveAllType<T>>::value, EigenArrayTag);
        DEFINE_DATA_CATEGORY(IsEigenMatrix<RemoveAllType<T>>::value, EigenMatrixTag);

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
        struct DataTraits<T, EigenArrayTag> : public DataTraits<T, EigenDenseTag>{};

        template <class T>
        struct DataTraits<T, EigenMatrixTag> : public DataTraits<T, EigenDenseTag>{};


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

            // eval
            template <class D>
            std::decay_t<typename Eigen::DenseBase<D>::EvalReturnType> Eval(const Eigen::DenseBase<D> & d){
                return d.eval();
            }

            // fill with scalar
            template <class D>
            std::decay_t<typename Eigen::DenseBase<D>::EvalReturnType> FillWithScalar(const Eigen::DenseBase<D> & d, typename D::Scalar s){
                auto de = d.eval();
                de.setConstant(s);
                return de;
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


            // cwise prod
            template <class D1, class D2>
            auto CWiseProd(const Eigen::MatrixBase<D1> & a, const Eigen::MatrixBase<D2> & b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b);
            }

            template <class D1, class D2>
            auto CWiseProd(const Eigen::ArrayBase<D1> & a, const Eigen::ArrayBase<D2> & b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b);
            }

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


            // sum elements
            template <class D>
            typename D::Scalar SumElements(const Eigen::DenseBase<D> & d) {
                return d.sum();
            }

        }


        namespace specific {

            // inverse matrix
            template <class D>
            auto InverseMatrix(const Eigen::MatrixBase<D> & m) -> decltype(m.inverse()) {
                return m.inverse();
            }        

        }



    }
}
 
#endif