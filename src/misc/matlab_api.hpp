#pragma once

#include "../core/basic_types.hpp"

namespace pano {
    namespace misc {

        class Matlab;
        class MAT;

        class MXA {
        public:
            MXA();
            MXA(void * mxa, bool dos = false);
            MXA(MXA && a);
            MXA & operator = (MXA && a);
            MXA(const MXA &) = delete;
            MXA & operator = (const MXA & a) = delete;
            virtual ~MXA();

        public:
            MXA clone(bool dos = false) const;
            void * mxa() const { return _mxa; }
            void * data() const;
            void setData(void* d);

#define DECL_MXARRAY_MEMBERFUNCTION_IS(what) bool is##what() const;

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

            size_t m() const;
            size_t n() const;
            void setM(size_t m);
            void setN(size_t n);

            bool null() const { return _mxa == 0; }
            bool operator! () const { return null(); }
            bool empty() const;

            bool isFromGlobalWorkspace() const;
            void setIsFromGlobalWorkspace(bool b);


            MXA(const cv::Mat & m, bool dos = false);
            MXA(double scalar, bool dos = false);
            MXA(const std::string & string, bool dos = false);
            MXA(const cv::SparseMat & m, bool dos = false);
            MXA(cv::InputArray m, bool dos = false);


            double scalar() const;
            std::string toString() const;
                       
            bool toCVOutputArray(cv::OutputArray o, bool lastDimIsChannel = true) const;
            cv::Mat toCVMat(bool lastDimIsChannel = true) const {
                cv::Mat m; toCVOutputArray(m, lastDimIsChannel); return m;
            }

            operator float() const { return (float)scalar(); }
            operator double() const { return scalar(); }
            operator long double() const { return scalar(); }

            operator int8_t() const { return static_cast<int8_t>(std::round(scalar())); }
            operator int16_t() const { return static_cast<int16_t>(std::round(scalar())); }
            operator int32_t() const { return static_cast<int32_t>(std::round(scalar())); }
            operator int64_t() const { return static_cast<int64_t>(std::round(scalar())); }

            operator uint8_t() const { return static_cast<uint8_t>(std::round(scalar())); }
            operator uint16_t() const { return static_cast<uint16_t>(std::round(scalar())); }
            operator uint32_t() const { return static_cast<uint32_t>(std::round(scalar())); }
            operator uint64_t() const { return static_cast<uint64_t>(std::round(scalar())); }

            operator std::string() const { return toString(); }
            operator cv::Mat() const { return toCVMat(); }

            size_t nelements() const;
            size_t nzmax() const;

            MXA cell(size_t i) const;
            void setCell(size_t i, const MXA & a);

            int nfields() const;
            const char * fieldName(int n) const;
            int fieldNumber(const std::string & name) const;
            MXA field(const std::string & name, int i = 0) const;
            void setField(const std::string & name, int i, const MXA & a);

            MXA property(const std::string & name, size_t i) const;
            void setProperty(const std::string & name, int i, const MXA & a);

            const char * className() const;

            size_t ndims() const;
            std::vector<size_t> dims() const;
            size_t dim(int d) const { return dims().at(d); }
            size_t length() const;

            size_t calcSingleSubscript(size_t dim, size_t * subs) const;

            template <class ... Ints>
            size_t offset(Ints... subs) const {
                size_t sbs[sizeof...(Ints)] = { subs ... };
                return calcSingleSubscript(sizeof...(Ints), sbs);
            }

            template <class T = double, class ... Ints, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            const T & at(Ints... subs) const {
                return static_cast<const T*>(data())[offset(subs...)];
            }
            template <class T = double, class ... Ints, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            T & at(Ints... subs) {
                return static_cast<T*>(data())[offset(subs...)];
            }

        protected:
            void * _mxa;
            bool _destroyWhenOutofScope;
        };


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
            bool null() const { return _fp == 0; }

            std::vector<std::string> varNames() const;

            MXA var(const std::string & name) const;
            bool setVar(const std::string & name, const MXA & mxa, bool asGlobal = false);
            bool removeVar(const std::string & name);

        private:
            std::string _fname;
            void * _fp;
        };


        // the matlab engine
        class Matlab {
        public:
            Matlab(const std::string & defaultDir = std::string(), bool singleUse = false);
            ~Matlab();

            Matlab(Matlab && e);
            Matlab & operator = (Matlab && e);

            Matlab(const Matlab &) = delete;
            Matlab & operator = (const Matlab &) = delete;

        public:
            bool started() const;
            bool run(const std::string & cmd) const;
            std::string lastMessage() const;
            bool errorLastRun() const;

            MXA var(const std::string & name) const;
            bool setVar(const std::string & name, const MXA & mxa);
            
            const Matlab & operator << (const std::string & cmd) const;

            bool cdAndAddAllSubfolders(const std::string & dir);

        private:
            char * _buffer;
            void * _eng;
        };


    }
}
