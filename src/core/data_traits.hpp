#ifndef PANOPTIC_CORE_DATA_TRAITS_HPP
#define PANOPTIC_CORE_DATA_TRAITS_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {
        
        static const int DynamicSizeValue = -1;

        namespace std_traits {
            namespace {
                template <typename T> struct IsComplexNumber : public std::false_type {};
                template <typename T> struct IsComplexNumber<std::complex<T>> : public std::is_floating_point<T>{};
            }

            template <typename T> struct IsActive : public std::integral_constant<bool,
                std::is_floating_point<T>::value || IsComplexNumber<T>::value> 
            {};

            template <typename T> using ElementType = T;
            template <typename T> struct IsLikeMatrix : public std::false_type {};
            template <typename T> struct StaticShape {
                static const int Rows = 1;
                static const int Cols = 1;
                static const int Size = 1;
            };

            // shape values
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetRows(const T & t) { return 1; }
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetCols(const T & t) { return 1; }
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetSize(const T & t) { return 1; }

            // shape operation
            template <typename T> std::enable_if_t<IsActive<T>::value> SetRows(T & t, int rows) {}
            template <typename T> std::enable_if_t<IsActive<T>::value> SetCols(T & t, int cols) {}
            template <typename T> std::enable_if_t<IsActive<T>::value> Resize(T & t, int rows, int cols) {}

            // evaluated type
            template <typename T> using EvaluatedType = std::decay_t<T>;

            // other operations
            // cast
            template <typename T> std::enable_if_t<IsActive<T>::value, T> data_cast(T && d) { return d; }

            

            // fill with scalar
            template <typename T> std::enable_if_t<IsActive<T>::value> Fill(T & t, const ElementType<T> & s) { t = s; }

            // element operation
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsSum(const T & t, ElementType<T> & s) { s = t; }
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsProd(const T & t, ElementType<T> & s) { s = t; }
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsMean(const T & t, ElementType<T> & s) { s = t; }

            // transpose
            template <typename T> using TransposeResultType = T;
            template <typename T> std::enable_if_t<IsActive<T>::value> Transpose(const T & from, TransposeResultType<T> & to) { to = from; }

            // elementwise op
            template <typename T> std::enable_if_t<IsActive<T>::value> ElementWiseMult(const T & a, const T & b, T & c) { c = a * b; }
            
            // matrix mult
            template <typename T1, typename T2> using MatrixMultResultType = T1;
            template <typename T1, typename T2> std::enable_if_t<IsActive<T1>::value && IsActive<T2>::value> 
                MatrixMult(const T1 & a, const T2 & b, MatrixMultResultType<T1, T2> & c) { c = a * b; }
        }



        namespace eigen_traits {

            using namespace Eigen;

            // is an eigen specified type, common types ingored
            template <typename T> struct IsActive : public std::integral_constant<bool, 
                std::is_base_of<DenseBase<std::decay_t<T>>, std::decay_t<T>>::value || 
                std::is_base_of<SparseMatrixBase<std::decay_t<T>>, std::decay_t<T>>::value>
            {};

            // evaluated type
            enum class EvaluatedTag { Dense, SparseMatrix, None };
            template <typename T, EvaluatedTag Flag = EvaluatedTag::None> struct Evaluated { using Type = std::decay_t<T>; };
            template <typename T> struct Evaluated<T, EvaluatedTag::Dense> { using Type = std::decay_t<decltype(std::declval<T>().eval())>; };
            template <typename T> struct Evaluated<T, EvaluatedTag::SparseMatrix> { using Type = std::decay_t<decltype(std::declval<T>().eval())>; };
            template <typename T> using EvaluatedType = 
                typename Evaluated<T, (
                std::is_base_of<DenseBase<std::decay_t<T>>, std::decay_t<T>>::value ? EvaluatedTag::Dense :
                std::is_base_of<SparseMatrixBase<std::decay_t<T>>, std::decay_t<T>>::value ? EvaluatedTag::SparseMatrix : 
                EvaluatedTag::None ) >::Type;

            // get element type
            template <typename T> struct Element { using Type = T; };
            template <typename T, int Rows, int Cols> struct Element<Matrix<T, Rows, Cols>>{ using Type = T; };
            template <typename T> struct Element<SparseMatrix<T>>{ using Type = T; };
            template <typename T> using ElementType = typename Element<T>::Type;

            // default multiplication is like matrix
            template <typename T> struct IsLikeMatrix : std::false_type {};
            template <typename T, int Rows, int Cols> struct IsLikeMatrix<Matrix<T, Rows, Cols>> : public std::true_type{};
            template <typename T> struct IsLikeMatrix<SparseMatrix<T>> : public std::true_type{};
            
            // static shape
            template <typename T> struct StaticShape {
                static const int Rows = 1;
                static const int Cols = 1;
                static const int Size = 1;
            };
            template <typename T, int R, int C> struct StaticShape<Matrix<T, R, C>> {
                static const int Rows = R == Eigen::Dynamic ? DynamicSizeValue : R;
                static const int Cols = C == Eigen::Dynamic ? DynamicSizeValue : C;
                static const int Size = (R == Eigen::Dynamic || C == Eigen::Dynamic) ? DynamicSizeValue : R * C;
            };
            template <typename T> struct StaticShape<SparseMatrix<T>> {
                static const int Rows = DynamicSizeValue;
                static const int Cols = DynamicSizeValue;
                static const int Size = DynamicSizeValue;
            };

            // shape values
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetRows(const T & t) { return t.rows(); }
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetCols(const T & t) { return t.cols(); }
            template <typename T> std::enable_if_t<IsActive<T>::value, int> GetSize(const T & t) { return t.size(); }

            // shape operation
            template <typename T> std::enable_if_t<IsActive<T>::value> SetRows(T & t, int rows) { t.resize(rows, t.cols()); }
            template <typename T> std::enable_if_t<IsActive<T>::value> SetCols(T & t, int cols) { t.resize(t.rows(), cols); }
            template <typename T> std::enable_if_t<IsActive<T>::value> Resize(T & t, int rows, int cols) { t.resize(rows, cols); }


            // other operations
            // fill with scalar
            template <typename Derived> void Fill(DenseBase<Derived> & t, const ElementType<Derived> & s) { t.fill(s); }
            template <typename Derived> void Fill(SparseMatrixBase<Derived> & t, const ElementType<Derived> & s) {
                if (s == 0){ t.derived().setZero(); }
                else if (t.size() == 0) {}
                else if (t.size() == 1) {
                    using TripletType = Eigen::Triplet<ElementType<Derived>, int>;
                    std::vector<TripletType> ts;
                    ts.reserve(t.size());
                    for (int i = 0; i < t.cols(); i++){
                        for (int j = 0; j < t.rows(); j++)
                            ts.push_back(TripletType(j, i, s));
                    }
                    t.derived().setFromTriplets(std::begin(ts), std::end(ts));
                }
            }

            // element operations
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsSum(const T & t, ElementType<T> & s) { s = t.sum(); }
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsProd(const T & t, ElementType<T> & s) { s = t.prod(); }
            template <typename T> std::enable_if_t<IsActive<T>::value> GetElementsMean(const T & t, ElementType<T> & s) { s = t.mean(); }

            // transpose
            template <typename T> struct TransposeResult { using Type = T; };
            template <typename T, int Rows, int Cols> struct TransposeResult<Matrix<T, Rows, Cols>>{ using Type = Matrix<T, Cols, Rows>; };
            template <typename T> using TransposeResultType = typename TransposeResult<T>::Type;

            template <typename T> std::enable_if_t<IsActive<T>::value> Transpose(const T & from, TransposeResultType<T> & to) {
                if (intptr_t(std::addressof(to)) == intptr_t(std::addressof(from))){ to.transposeInPlace(); } 
                else { to = from.transpose(); }
            }

            // elementwise mult
            template <typename T> std::enable_if_t<IsActive<T>::value> ElementWiseMult(const T & a, const T & b, T & c) { c = a.cwiseProduct(b); }
            
            // matrix mult
            template <typename T1, typename T2> struct MatrixMultResult {};
            template <typename T, int A, int B, int C> struct MatrixMultResult<Matrix<T, A, B>, Matrix<T, B, C>> { using Type = Matrix<T, A, C>; };
            template <typename T> struct MatrixMultResult<SparseMatrix<T>, SparseMatrix<T>> { using Type = SparseMatrix<T>; };
            template <typename T1, typename T2> using MatrixMultResultType = typename MatrixMultResult<T1, T2>::Type;
            template <typename T1, typename T2> std::enable_if_t<IsActive<T1>::value && IsActive<T2>::value>
                MatrixMult(const T1 & a, const T2 & b, MatrixMultResultType<T1, T2> & ab) { ab = a * b; }

            // exp
            template <typename T> std::enable_if_t<IsActive<T>::value, T> exp(const T & t) { return t.array().exp(); }
        }
        
        
        

        template <typename T> struct IsActive : public std::integral_constant<bool, 
            std_traits::IsActive<T>::value || 
            eigen_traits::IsActive<T>::value > 
        {};

        template <typename T> using ElementType = 
            std::conditional_t <std_traits::IsActive<T>::value, std_traits::ElementType<T>, 
            std::conditional_t <eigen_traits::IsActive<T>::value, eigen_traits::ElementType<T>, 
            void>>;

        template <typename T> struct IsLikeMatrix : public
            std::conditional_t <std_traits::IsActive<T>::value, std_traits::IsLikeMatrix<T>,
            std::conditional_t <eigen_traits::IsActive<T>::value, eigen_traits::IsLikeMatrix<T>,
            std::false_type>>
        {};

        template <typename T> struct StaticShape : public
            std::conditional_t <std_traits::IsActive<T>::value, std_traits::StaticShape<T>,
            std::conditional_t <eigen_traits::IsActive<T>::value, eigen_traits::StaticShape<T>,
            std::false_type>>
        {};


        template <typename T> using EvaluatedType = 
            std::conditional_t <eigen_traits::IsActive<T>::value, eigen_traits::EvaluatedType<T>,
            std::conditional_t <std_traits::IsActive<T>::value, std_traits::EvaluatedType<T>,
            void >> ;

        
        using std_traits::GetRows;
        using std_traits::GetCols;
        using std_traits::GetSize;

        using eigen_traits::GetRows;
        using eigen_traits::GetCols;
        using eigen_traits::GetSize;

        using std_traits::SetRows;
        using std_traits::SetCols;
        using std_traits::Resize;

        using eigen_traits::SetRows;
        using eigen_traits::SetCols;
        using eigen_traits::Resize;


        using std_traits::Fill;
        using eigen_traits::Fill;

        using std_traits::GetElementsSum;
        using std_traits::GetElementsProd;
        using std_traits::GetElementsMean;

        using eigen_traits::GetElementsSum;
        using eigen_traits::GetElementsProd;
        using eigen_traits::GetElementsMean;

        template <typename T> using TransposeResultType =
            std::conditional_t <std_traits::IsActive<T>::value, std_traits::TransposeResultType<T>,
            std::conditional_t <eigen_traits::IsActive<T>::value, eigen_traits::TransposeResultType<T>,
            void>>;

        using std_traits::Transpose;
        using eigen_traits::Transpose;

        using std_traits::ElementWiseMult;
        using eigen_traits::ElementWiseMult;

        template <typename T1, typename T2> using MatrixMultResultType =
            std::conditional_t <std_traits::IsActive<T1>::value && std_traits::IsActive<T2>::value, 
            std_traits::MatrixMultResultType<T1, T2>,
            std::conditional_t <eigen_traits::IsActive<T1>::value && eigen_traits::IsActive<T2>::value, 
            eigen_traits::MatrixMultResultType<T1, T2>,
            void >> ;

        using std_traits::MatrixMult;
        using eigen_traits::MatrixMult;

        template <typename T> struct IsStaticallyScalar : public std::integral_constant<bool, 
            StaticShape<T>::Size == 1> 
        {};

        using std::exp;
        using eigen_traits::exp;

        template <typename To, typename From> const To & ValueOf(const From & f) { return f; }
        template <typename T> std::enable_if_t<std::is_pointer<T>::value, const std::remove_pointer_t<T> &> 
            ValueOf(T const & f) { return *f };


    }
}
 
#endif