#ifndef PANORAMIX_SANDBOX_FOR_EXPRESSION2_HPP
#define PANORAMIX_SANDBOX_FOR_EXPRESSION2_HPP

#include <Eigen/Dense>

#include "for_expression.hpp"

namespace panoramix {
    namespace sandbox {

        DEFINE_CATEGORY(std::is_floating_point<T>::value, floating_point_tag);
        template <class T>
        struct Traits<T, floating_point_tag> {
            static const bool c = true;
            static void print() {
                std::cout << "is floating point!" << std::endl;
            }
        };

        DEFINE_CATEGORY(std::is_integral<T>::value, integral_tag);
        template <class T>
        struct Traits<T, integral_tag> {
            static void print() {
                std::cout << "is integer!" << std::endl;
            }
        };

        namespace {

            template <class T>
            struct is_eigen_dense
                : public std::is_base_of<Eigen::EigenBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_dense<Eigen::DenseBase<Derived>> : public std::true_type{};

            template <class T>
            struct is_eigen_matrix
                : public std::is_base_of<Eigen::MatrixBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_matrix<Eigen::MatrixBase<Derived>> : public std::true_type{};

            template <class T>
            struct is_eigen_array
                : public std::is_base_of<Eigen::ArrayBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_array<Eigen::ArrayBase<Derived>> : public std::true_type{};

        }

        DEFINE_CATEGORY(is_eigen_dense<T>::value && !is_eigen_array<T>::value && !is_eigen_matrix<T>::value, eigen_dense_tag);
        template <class T>
        struct Traits<T, eigen_dense_tag> {
            static void print() {
                std::cout << "is eigen dense!" << std::endl;
            }
        };


        DEFINE_CATEGORY(is_eigen_matrix<T>::value, eigen_matrix_tag);
        template <class T>
        struct Traits<T, eigen_matrix_tag> {
            static void print() {
                std::cout << "is eigen matrix!" << std::endl;
            }
        };

        DEFINE_CATEGORY(is_eigen_array<T>::value, eigen_array_tag);
        template <class T>
        struct Traits<T, eigen_array_tag> {
            static void print() {
                std::cout << "is eigen array!" << std::endl;
            }
        };

    }
}
 
#endif