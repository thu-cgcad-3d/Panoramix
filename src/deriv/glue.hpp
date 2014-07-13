#ifndef PANORAMIX_DERIV_GLUE_HPP
#define PANORAMIX_DERIV_GLUE_HPP

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
 
namespace panoramix {
    namespace deriv {

        template <class T, int M, int N>
        Eigen::Matrix<T, M, N> CVMatToEigenMat(const cv::Matx<T, M, N> & m) {
            Eigen::Matrix<T, M, N> mm;
            for (int i = 0; i < M; i++){
                for (int j = 0; j < N; j++){
                    mm(i, j) = m(i, j);
                }
            }
            return mm;
        }

        template <class T, int M, int N>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> CVMatToEigenMatX(const cv::Matx<T, M, N> & m) {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mm;
            for (int i = 0; i < M; i++){
                for (int j = 0; j < N; j++){
                    mm(i, j) = m(i, j);
                }
            }
            return mm;
        }

        template <class T, int M, int N, int O, int MM, int NN>
        cv::Matx<T, M, N> EigenMatToCVMat(const Eigen::Matrix<T, M, N, O, MM, NN> & m) {
            cv::Matx<T, M, N> mm;
            for (int i = 0; i < M; i++){
                for (int j = 0; j < N; j++){
                    mm(i, j) = m(i, j);
                }
            }
            return mm;
        }

        template <class T, int M, int O, int MM, int NN>
        cv::Vec<T, M> EigenVecToCVVec(const Eigen::Matrix<T, M, 1, O, MM, NN> & m) {
            cv::Vec<T, M> mm;
            for (int i = 0; i < M; i++)
                mm(i) = m(i);
            return mm;
        }

    }
}
 
#endif