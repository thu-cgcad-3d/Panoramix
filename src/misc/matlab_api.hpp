#ifndef PANORAMIX_MISC_MXARRAY_HPP
#define PANORAMIX_MISC_MXARRAY_HPP

#ifndef USE_MATLAB
#error USE_MATLAB flag required
#endif

#include <mex.h>
#include <mat.h>

#include "../core/basic_types.hpp"

namespace panoramix {
    namespace misc {

        class MXArray {
        public:
            MXArray();
            MXArray(mxArray * mxa, bool dos = false);
            MXArray(MXArray && a);
            MXArray & operator = (MXArray && a);
            MXArray(const MXArray &) = delete;
            MXArray & operator = (const MXArray & a) = delete;
            virtual ~MXArray();

        public:
            operator mxArray * () const { return _mxa; }
            mxClassID classID() const { return mxGetClassID(_mxa); }
            void * data() const { return mxGetData(_mxa); }
            void setData(void* d) { mxSetData(_mxa, d); }

#define DECL_MXARRAY_MEMBERFUNCTION_IS(what) bool is##what() const { return mxIs##what(_mxa); }

            DECL_MXARRAY_MEMBERFUNCTION_IS(Numeric)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Cell)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Logical)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Char)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Struct)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Opaque)
            DECL_MXARRAY_MEMBERFUNCTION_IS(FunctionHandle)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Object)

            DECL_MXARRAY_MEMBERFUNCTION_IS(Complex)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Sparse)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Double)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Single)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Int8)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Uint8)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Int16)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Uint16)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Int32)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Uint32)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Int64)
            DECL_MXARRAY_MEMBERFUNCTION_IS(Uint64)


            template <class T> bool is() const { return false; }
            template <> bool is<double>() const { return isDouble(); }
            template <> bool is<float>() const { return isSingle(); }
            template <> bool is<int8_t>() const { return isInt8(); }
            template <> bool is<uint8_t>() const { return isUint8(); }
            template <> bool is<int16_t>() const { return isInt16(); }
            template <> bool is<uint16_t>() const { return isUint16(); }
            template <> bool is<int32_t>() const { return isInt32(); }
            template <> bool is<uint32_t>() const { return isUint32(); }
            template <> bool is<int64_t>() const { return isInt64(); }
            template <> bool is<uint64_t>() const { return isUint64(); }

            size_t m() const { return mxGetM(_mxa); }
            size_t n() const { return mxGetN(_mxa); }
            void setM(size_t m) { mxSetM(_mxa, m); }
            void setN(size_t n) { mxSetN(_mxa, n); }

            bool null() const { return _mxa == nullptr; }
            bool operator! () const { return null(); }
            operator bool() const { return !null(); }
            bool empty() const { return mxIsEmpty(_mxa); }

            bool isFromGlobalWorkspace() const { return mxIsFromGlobalWS(_mxa); }
            void setIsFromGlobalWorkspace(bool b) { mxSetFromGlobalWS(_mxa, b); }

            double scalar() const { return mxGetScalar(_mxa); }
            std::string toString() const {
                size_t len = nelements() + 1;
                char * buffer = new char[len];
                memset(buffer, '\0', len);
                mxGetString(_mxa, buffer, len);
                std::string str = buffer;
                delete[] buffer;
                return str;
            }

            operator double() const { return scalar(); }
            operator std::string() const { return toString(); }

            size_t nelements() const { return mxGetNumberOfElements(_mxa); }
            size_t nzmax() const { return mxGetNzmax(_mxa); }

            MXArray cell(size_t i) const { return mxGetCell(_mxa, i); }
            void setCell(size_t i, const MXArray & a) { mxSetCell(_mxa, i, a._mxa); }

            int nfields() const { return mxGetNumberOfFields(_mxa); }
            const char * fieldName(int n) const { return mxGetFieldNameByNumber(_mxa, n); }
            int fieldNumber(const std::string & name) const { return mxGetFieldNumber(_mxa, name.c_str()); }
            MXArray field(const std::string & name, int i) const { return mxGetField(_mxa, i, name.c_str()); }
            void setField(const std::string & name, int i, const MXArray & a) { mxSetField(_mxa, i, name.c_str(), a._mxa); }

            MXArray property(const std::string & name, size_t i) const { return mxGetProperty(_mxa, i, name.c_str()); }
            void setProperty(const std::string & name, int i, const MXArray & a) { mxSetProperty(_mxa, i, name.c_str(), a._mxa); }

            const char * className() const { return mxGetClassName(_mxa); }

            size_t ndims() const { return mxGetNumberOfDimensions(_mxa); }
            std::vector<size_t> dims() const {
                std::vector<size_t> ds(ndims());
                std::copy_n(mxGetDimensions(_mxa), ds.size(), ds.begin());
                return ds;
            }
            size_t dim(int d) const {
                return dims().at(d);
            }
            size_t length() const {
                auto ds = dims();
                if (ds.empty()) return 0;
                return *std::max_element(ds.begin(), ds.end());
            }

            template <class ... Ints>
            size_t offset(Ints... subs) const {
                size_t sbs[sizeof...(Ints)] = { subs ... };
                return mxCalcSingleSubscript(_mxa, sizeof...(Ints), sbs);
            }

            template <class T, class ... Ints, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            const T & at(Ints... subs) const {
                return static_cast<const T*>(mxGetData(_mxa))[offset(subs...)];
            }
            template <class T, class ... Ints, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            T & at(Ints... subs) {
                return static_cast<T*>(mxGetData(_mxa))[offset(subs...)];
            }


        protected:
            mxArray * _mxa;
            bool _destroyWhenOutofScope;
        };


        inline MXArray ToMXArray(mxArray * mxa, bool dos = false) {
            return MXArray(mxa, dos);
        }

        inline const MXArray ToMXArray(const mxArray * mxa){
            return MXArray(const_cast<mxArray*>(mxa), false);
        }

        inline std::vector<MXArray> ToMXArray(int n, mxArray ** mxas, bool dos = false){
            std::vector<MXArray> as;
            for (int i = 0; i < n; i++){
                as.emplace_back(mxas[i], dos);
            }
            return as;
        }

        inline std::vector<const MXArray> ToMXArray(int n, const mxArray ** mxas){
            std::vector<const MXArray> as;
            for (int i = 0; i < n; i++){
                as.emplace_back(const_cast<mxArray*>(mxas[i]), false);
            }
            return as;
        }


        MXArray CVInputArrayToMXArray(cv::InputArray mat, bool dos = false);
        bool CVOutputArrayFromMXArray(const MXArray & mxa, cv::OutputArray mat, bool lastDimIsChannel = true);

        inline MXArray ToMXArray(const core::Image & mat, bool dos = false) {
            return CVInputArrayToMXArray(mat, dos); 
        }
        inline bool FromMXArray(const MXArray & mxa, core::Image & mat, bool lastDimIsChannel = true){ 
            return CVOutputArrayFromMXArray(mxa, mat, lastDimIsChannel); 
        }

        MXArray ToMXArray(const std::vector<cv::KeyPoint> & ps, bool dos = false);
        inline bool FromMXArray(const MXArray & mxa, std::vector<cv::KeyPoint> & ps) {
            return CVOutputArrayFromMXArray(mxa, ps, false); 
        }

        inline MXArray ToMXArray(double d, bool dos = false){ return CVInputArrayToMXArray(d, dos); }
        inline bool FromMXArray(const MXArray & mxa, double & d) {
            cv::Mat m;
            bool r = CVOutputArrayFromMXArray(mxa, m, false); 
            d = m.at<double>(0);
            return r;
        }

        MXArray ToMXArray(const cv::SparseMat & mat, bool dos = false);


        // the mat file
        class MAT {
        public:
            enum OpeningMode {
                Read,
                Update,
                Write,
                Write_4,
                Write_CompressedData,
                Write_7_3
            };

        public:
            MAT();
            explicit MAT(const std::string & fname, const std::string & mode);
            explicit MAT(const std::string & fname, OpeningMode mode);
            MAT(MAT && a);
            MAT & operator = (MAT && a);
            MAT(const MAT &) = delete;
            MAT & operator = (const MAT & a) = delete;
            virtual ~MAT();

        public:
            bool null() const { return _fp == nullptr; }

            std::vector<std::string> varNames() const;

            MXArray var(const std::string & name) const { 
                return matGetVariable(_fp, name.c_str());
            }
            bool setVar(const std::string & name, const MXArray & mxa, bool asGlobal = false) { 
                if (!asGlobal){
                    return matPutVariable(_fp, name.c_str(), mxa) == 0;
                }
                else{
                    return matPutVariableAsGlobal(_fp, name.c_str(), mxa) == 0;
                }
            }
            bool removeVar(const std::string & name) {
                return matDeleteVariable(_fp, name.c_str()) == 0;
            }

        private:
            std::string _fname;
            MATFile * _fp;
        };

    }
}

#endif