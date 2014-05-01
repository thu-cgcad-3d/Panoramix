#ifndef PANORAMIX_CORE_DATA_TRAITS_HPP
#define PANORAMIX_CORE_DATA_TRAITS_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <armadillo>

#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        namespace {

            template <class T>
            struct is_eigen_dense 
                : public std::is_base_of<Eigen::EigenBase<T>, T>
            {};

            template <class T>
            struct is_eigen_matrix
                : public std::is_base_of<Eigen::MatrixBase<T>, T>
            {};

            template <class T>
            struct is_eigen_array
                : public std::is_base_of<Eigen::ArrayBase<T>, T>
            {};


           /* template <class T>
            struct is_arma_mat<arma::Mat<T>> : public std::true_type {};*/

        }

        struct floating_point_tag {};
        struct eigen_dense_tag {};
        struct eigen_matrix_tag {};
        struct eigen_array_tag {};
        struct arma_mat_tag {};

        template <class T>
        struct type_tag {
            using type = std::conditional_t < std::is_floating_point<std::decay_t<T>>::value, floating_point_tag,
                std::conditional_t < is_eigen_matrix<std::decay_t<T>>::value, eigen_matrix_tag,
                std::conditional_t < is_eigen_array<std::decay_t<T>>::value, eigen_array_tag,
                std::conditional_t < is_eigen_dense<std::decay_t<T>>::value, eigen_dense_tag,
                //std::conditional_t < is_arma_mat<std::decay_t<T>>::value, arma_mat_tag,
                void >> >> ;
        };
                
        template <class T>
        using type_tag_t = typename type_tag<T>::type;

        // type is currently not supported
        template <class T, typename Tag = typename type_tag<T>::type >
        struct is_not_supported : public std::false_type {};
        template <class T>
        struct is_not_supported<T, void> : public std::true_type{};


        // type traits
        template <class T, class Tag = typename type_tag<T>::type>
        struct traits {};
        

        // floating points
        template <class T>
        struct traits<T, floating_point_tag> {

            static const bool is_scalar = true;

            // the type for storage
            using storage_type = std::decay_t<T>;

            // scalar type
            using scalar_type = storage_type;

            // fill
            static storage_type fill(const storage_type & t, scalar_type s) { return storage_type(s); }

            // should be cached
            static const bool should_be_cached = false;

            // eval
            static storage_type eval(const storage_type & t) { return t; }

            // get one
            static storage_type one() { return static_cast<storage_type>(1.0); }

            // get zero
            static storage_type zero() { return static_cast<storage_type>(0); }

            // transpose
            static T transpose(T t) { return t; }

            // pow
            static storage_type pow(T t, scalar_type exponent) { return std::pow(t, exponent); }

        };

        // eigen_dense
        template <class T>
        struct traits<T, eigen_dense_tag> {

            static const bool is_scalar = false;

            // the type for storage
            using storage_type = std::decay_t<decltype(std::declval<T>().eval())>;

            // scalar type
            using scalar_type = typename storage_type::Scalar;

            // fill
            static storage_type fill(const storage_type & t, scalar_type s) { 
                //static_assert(false, "1");
                storage_type tt = t;
                tt.setConstant(s);
                return tt;
            }

            // should be cached
            static const bool should_be_cached = std::is_same<T, storage_type>::value;
            //static const bool should_be_cached = false;

            // eval
            static storage_type eval(const storage_type & t) { return t.eval(); }

            // get one
            static storage_type one() { assert(false && "cannot get one from eigen dense (which is not a scalar) type!"); throw ""; }
            
            // get zero
            // this is not implemented since no shape information is given for non scalar data (e.g MatrixXd)
            static storage_type zero() { assert(false && "cannot get zero from eigen dense (which is not a scalar) type!"); throw ""; }

        };

        // eigen matrix
        template <class T>
        struct traits<T, eigen_matrix_tag> : public traits<T, eigen_dense_tag> { 

            // transpose
            static auto transpose(T t) -> decltype(t.transpose()) { return t.transpose(); }

            // pow
            static auto pow(T t, scalar_type exponent) -> decltype(t.array().pow(exponent).matrix().eval()) { 
                return t.array().pow(exponent).matrix().eval();
            }

        };

        // eigen array
        template <class T>
        struct traits<T, eigen_array_tag> : public traits<T, eigen_dense_tag>{ 

            // transpose
            static T transpose(T t) { return t; }

            // pow
            static auto pow(T t, scalar_type exponent) -> decltype(t.pow(exponent)) {
                return t.pow(exponent); 
            }

        };

        // alias types
        template <class T>
        using storage_type = typename traits<T>::storage_type;
        template <class T>
        struct is_storage_type : public std::is_same<T, storage_type<T>> {};
        template <class T>
        using scalar_type = typename traits<T>::scalar_type;
        template <class T>
        struct is_scalar_type : public std::is_same<T, scalar_type<T>> {};




        // type pair traits
        template <class T1, class T2, 
        class Tag1 = typename type_tag<T1>::type, 
        class Tag2 = typename type_tag<T2>::type>
        struct traits2 {};

        template <class T1, class T2>
        struct traits2<T1, T2, floating_point_tag, floating_point_tag> {            
            static T1 cast(T2 from)  { return static_cast<T1>(from); } 
            static auto cwise_product(T1 && a, T2 && b) -> decltype(a*b) {
                return a * b;
            }
        };

        template <class T1, class T2>
        struct traits2<T1, T2, eigen_dense_tag, eigen_dense_tag> {
        };


        template <class T1, class T2>
        struct traits2<T1, T2, eigen_matrix_tag, eigen_matrix_tag> 
            : public traits2<T1, T2, eigen_dense_tag, eigen_dense_tag> {
            static T1 cast(T2 from)  { 
                return from.cast<typename T1::Scalar>().eval(); 
            }
            static auto cwise_product(T1 && a, T2 && b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b); 
            }
        };

        template <class T1, class T2>
        struct traits2<T1, T2, eigen_array_tag, eigen_array_tag> 
            : public traits2<T1, T2, eigen_dense_tag, eigen_dense_tag>{
            static T1 cast(T2 from)  {
                return from.cast<typename T1::Scalar>().eval();
            }
            static auto cwise_product(T1 && a, T2 && b) -> decltype(a.cwiseProduct(b)) {
                return a.cwiseProduct(b);
            }
        };

        template <class T1, class T2>
        struct traits2<T1, T2, eigen_array_tag, eigen_matrix_tag> 
            : public traits2<T1, T2, eigen_dense_tag, eigen_dense_tag>{
            static T1 cast(T2 from) {
                return from.array().cast<typename T1::Scalar>().eval();
            }
            static auto cwise_product(T1 && a, T2 && b) -> decltype(a.cwiseProduct(b.array())) {
                return a.cwiseProduct(b.array());
            }
        };

        template <class T1, class T2>
        struct traits2<T1, T2, eigen_matrix_tag, eigen_array_tag> 
            : public traits2<T1, T2, eigen_dense_tag, eigen_dense_tag>{
            static T1 cast(T2 from) {
                return from.matrix().cast<typename T1::Scalar>().eval();
            }
            static auto cwise_product(T1 && a, T2 && b) -> decltype(a.cwiseProduct(b.matrix())) {
                return a.cwiseProduct(b.matrix());
            }
        };


        // data cast
        template <class T1, class T2>
        inline T1 data_cast(T2 && from) {
            return traits2<T1, T2>::cast(from);
        }

        // cwise product
        template <class T1, class T2>
        inline auto data_cwise_product(T1 && a, T2 && b) 
            -> decltype(traits2<T1, T2>::cwise_product(a, b)) {
            static_assert(std::is_same<scalar_type<T1>, scalar_type<T2>>::value, "scalar types MUST be the same!");
            return traits2<T1, T2>::cwise_product(a, b);
        }

    }
}
 
#endif