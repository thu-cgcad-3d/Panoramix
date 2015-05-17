#include "matlab_api.hpp"

#ifndef USE_MATLAB
#error USE_MATLAB flag required
#endif

namespace panoramix {

    namespace misc {

        // template helpers
        template <typename T>	struct RMXType	{ static const mxClassID ClassID = mxUNKNOWN_CLASS; };
        template <int ID>		struct RMXID	{ using type = void; };

#define CONNECT_CLASSID(T, id) \
    template <> struct RMXType<T>	{ static const mxClassID ClassID = id; }; \
    template <> struct RMXID<id>	{ using type = T; }

        CONNECT_CLASSID(mxLogical, mxLOGICAL_CLASS);
        CONNECT_CLASSID(double, mxDOUBLE_CLASS);
        CONNECT_CLASSID(float, mxSINGLE_CLASS);
        CONNECT_CLASSID(int8_t, mxINT8_CLASS);
        CONNECT_CLASSID(int16_t, mxINT16_CLASS);
        CONNECT_CLASSID(int32_t, mxINT32_CLASS);
        CONNECT_CLASSID(int64_t, mxINT64_CLASS);
        CONNECT_CLASSID(uint8_t, mxUINT8_CLASS);
        CONNECT_CLASSID(uint16_t, mxUINT16_CLASS);
        CONNECT_CLASSID(uint32_t, mxUINT32_CLASS);
        CONNECT_CLASSID(uint64_t, mxUINT64_CLASS);

        template <typename T>
        inline mxArray* CreateNumericMatrix(const T* d, int m, int n) {
            mxArray* ma = mxCreateNumericMatrix(m, n, RMXType<T>::ClassID, mxREAL);
            std::memcpy(mxGetData(ma), d, sizeof(T)*m*n);
            return ma;
        }
        template <typename T>
        inline mxArray* CreateNumericMatrix(int m, int n){
            return mxCreateNumericMatrix(m, n, RMXType<T>::ClassID, mxREAL);
        }

        inline mxClassID CVDepthToMxClassID(int cvDepth){
            switch (cvDepth){
            case CV_8U: return mxUINT8_CLASS;
            case CV_8S: return mxINT8_CLASS;
            case CV_16U: return mxUINT16_CLASS;
            case CV_16S: return mxINT16_CLASS;
            case CV_32S: return mxINT32_CLASS;
            case CV_32F: return mxSINGLE_CLASS;
            case CV_64F: return mxDOUBLE_CLASS;
            default:
                std::cout << "this cv depth type cannot be converted to a matlab class id" << std::endl;
                return mxUNKNOWN_CLASS;
            }
        }

        inline int MxClassIDToCVDepth(mxClassID id){
            switch (id){
            case mxLOGICAL_CLASS: return CV_8U;
            case mxDOUBLE_CLASS: return CV_64F;
            case mxSINGLE_CLASS: return CV_32F;
            case mxINT8_CLASS: return CV_8S;
            case mxINT16_CLASS: return CV_16S;
            case mxINT32_CLASS: return CV_32S;
            case mxUINT8_CLASS: return CV_8U;
            case mxUINT16_CLASS: return CV_16U;
            case mxUINT32_CLASS: return CV_32S;
            default:
                std::cout << "this matlab class id cannot be converted to a cv depth type" << std::endl;
                return -1;
            }
        }


        MXArray::MXArray() : _mxa(nullptr), _destroyWhenOutofScope(false) {}
        MXArray::MXArray(mxArray * mxa, bool dos) : _mxa(mxa), _destroyWhenOutofScope(dos){}
        MXArray::MXArray(MXArray && a) {
            _mxa = a._mxa;
            a._mxa = nullptr;
            _destroyWhenOutofScope = a._destroyWhenOutofScope;
            a._destroyWhenOutofScope = false;
        }
        MXArray & MXArray::operator = (MXArray && a) {
            std::swap(_mxa, a._mxa);
            std::swap(_destroyWhenOutofScope, a._destroyWhenOutofScope);
            return *this;
        }
        MXArray::~MXArray() {
            if (_destroyWhenOutofScope){
                mxDestroyArray(_mxa);
            }
            _mxa = nullptr;
        }

        cv::Mat MXArray::toCVMat(bool lastDimIsChannel) const{
            if (_mxa == nullptr){
                return cv::Mat();
            }
            cv::Mat mat;
            const mxArray * ma = _mxa;
            int d = mxGetNumberOfDimensions(ma);
            assert((d >= 2) && "dimension num of the variable must be over 2");
            const mwSize * dimSizes = mxGetDimensions(ma);

            // set channels, dims and dimSizes
            int channels = 0;
            int cvDims = 0;
            int * cvDimSizes = nullptr;

            if (lastDimIsChannel && d > 2){
                channels = dimSizes[d - 1];
                cvDims = d - 1;
                cvDimSizes = new int[d - 1];
                std::copy(dimSizes, dimSizes + d - 1, cvDimSizes);
            }
            else{
                channels = 1;
                cvDims = d;
                cvDimSizes = new int[d];
                std::copy(dimSizes, dimSizes + d, cvDimSizes);
            }

            const size_t szForEachElem = mxGetElementSize(ma);
            int depth = MxClassIDToCVDepth(mxGetClassID(ma));
            if (depth == -1){
                delete[] cvDimSizes;
                return cv::Mat();
            }

            const uint8_t * mad = (const uint8_t*)mxGetData(ma);
            // create Mat
            mat.create(cvDims, cvDimSizes, CV_MAKETYPE(depth, channels));
            cv::Mat im = mat;

            mwIndex * mxIndices = new mwIndex[im.dims + 1];
            int * cvIndices = new int[im.dims];

            cv::MatConstIterator iter(&im);
            int imTotal = im.total();
            for (int i = 0; i < imTotal; i++, ++iter){
                // get indices in cv::Mat
                iter.pos(cvIndices);
                // copy indices to mxIndices
                std::copy(cvIndices, cvIndices + im.dims, mxIndices);
                for (mwIndex k = 0; k < channels; k++){
                    uint8_t * toDataHead = (*iter) + k * szForEachElem;
                    mxIndices[im.dims] = k; // set the last indices
                    const uint8_t * fromDataHead = mad + mxCalcSingleSubscript(ma, im.dims + 1, mxIndices) * szForEachElem;
                    std::memcpy(toDataHead, fromDataHead, szForEachElem);
                }
            }

            delete[] cvDimSizes;
            delete[] mxIndices;
            delete[] cvIndices;
            return mat;
        }


        MXArray CVInputArrayToMXArray(cv::InputArray mat, bool dos){
            cv::Mat& im = mat.getMat();

            // collect all dimensions of im
            int channelNum = im.channels();
            mwSize * dimSizes = new mwSize[im.dims + 1];
            for (int i = 0; i < im.dims; i++)
                dimSizes[i] = im.size[i];
            dimSizes[im.dims] = channelNum;

            // create mxArray
            mxArray* ma = mxCreateNumericArray(im.dims + 1, dimSizes, CVDepthToMxClassID(im.depth()), mxREAL);
            delete[] dimSizes;

            if (!ma)
                return nullptr;

            uint8_t * mad = (uint8_t*)mxGetData(ma);
            const size_t szForEachElem = im.elemSize1();
            mwIndex * mxIndices = new mwIndex[im.dims + 1];
            int * cvIndices = new int[im.dims];

            cv::MatConstIterator iter(&im);
            int imTotal = im.total();
            for (int i = 0; i < imTotal; i++, ++iter){
                // get indices in cv::Mat
                iter.pos(cvIndices);
                // copy indices to mxIndices
                std::copy(cvIndices, cvIndices + im.dims, mxIndices);
                for (mwIndex k = 0; k < channelNum; k++){
                    const uint8_t * fromDataHead = (*iter) + k * szForEachElem;
                    mxIndices[im.dims] = k; // set the last indices
                    uint8_t * toDataHead = mad + mxCalcSingleSubscript(ma, im.dims + 1, mxIndices) * szForEachElem;
                    std::memcpy(toDataHead, fromDataHead, szForEachElem);
                }
            }

            delete[] mxIndices;
            delete[] cvIndices;
            return MXArray(ma, dos);
        }

        bool CVOutputArrayFromMXArray(const MXArray & mxa, cv::OutputArray mat, bool lastDimIsChannel){
            if (!mat.needed())
                return true;
            if (mxa.null())
                return false;
            const mxArray * ma = mxa;
            int d = mxGetNumberOfDimensions(ma);
            assert((d >= 2) && "dimension num of the variable must be over 2");
            const mwSize * dimSizes = mxGetDimensions(ma);

            // set channels, dims and dimSizes
            int channels = 0;
            int cvDims = 0;
            int * cvDimSizes = nullptr;

            if (lastDimIsChannel && d > 2){
                channels = dimSizes[d - 1];
                cvDims = d - 1;
                cvDimSizes = new int[d - 1];
                std::copy(dimSizes, dimSizes + d - 1, cvDimSizes);
            }
            else{
                channels = 1;
                cvDims = d;
                cvDimSizes = new int[d];
                std::copy(dimSizes, dimSizes + d, cvDimSizes);
            }

            const size_t szForEachElem = mxGetElementSize(ma);
            int depth = MxClassIDToCVDepth(mxGetClassID(ma));
            if (depth == -1){
                delete[] cvDimSizes;
                return false;
            }

            const uint8_t * mad = (const uint8_t*)mxGetData(ma);
            // create Mat
            mat.create(cvDims, cvDimSizes, CV_MAKETYPE(depth, channels));
            cv::Mat im = mat.getMat();

            mwIndex * mxIndices = new mwIndex[im.dims + 1];
            int * cvIndices = new int[im.dims];

            cv::MatConstIterator iter(&im);
            int imTotal = im.total();
            for (int i = 0; i < imTotal; i++, ++iter){
                // get indices in cv::Mat
                iter.pos(cvIndices);
                // copy indices to mxIndices
                std::copy(cvIndices, cvIndices + im.dims, mxIndices);
                for (mwIndex k = 0; k < channels; k++){
                    uint8_t * toDataHead = (*iter) + k * szForEachElem;
                    mxIndices[im.dims] = k; // set the last indices
                    const uint8_t * fromDataHead = mad + mxCalcSingleSubscript(ma, im.dims + 1, mxIndices) * szForEachElem;
                    std::memcpy(toDataHead, fromDataHead, szForEachElem);
                }
            }

            delete[] cvDimSizes;
            delete[] mxIndices;
            delete[] cvIndices;
            return true;
        }

        MXArray ToMXArray(const std::vector<cv::KeyPoint> & ps, bool dos){
            MXArray ma = mxCreateDoubleMatrix(2, ps.size(), mxREAL);
            if (ma.null())
                return nullptr;
            for (int i = 0; i < ps.size(); i++){
                ma.at<double>(0, i) = ps[i].pt.x;
                ma.at<double>(1, i) = ps[i].pt.y;
            }
            return ma;
        }


        MXArray ToMXArray(const cv::SparseMat & mat, bool dos){
            auto & im = mat;

            // collect all dimensions of im
            int channelNum = im.channels();
            assert(channelNum == 1);

            // create mxArray
            int nzc = mat.nzcount();
            mxArray* ma = mxCreateSparse(im.size(0), im.size(1), nzc, mxREAL);

            if (!ma)
                return nullptr;

            std::vector<std::tuple<int, int, double>> triplets; // col - row - val
            triplets.reserve(nzc);

            THERE_ARE_BOTTLENECKS_HERE(); // we should only iterate non zero elements
            for (auto iter = mat.begin(); iter != mat.end(); ++iter){
                int ii = iter.node()->idx[0];
                int jj = iter.node()->idx[1];
                uchar* data = iter.ptr;
                double v = 0.0;
                if (mat.type() == CV_32FC1){
                    float value = 0.0f;
                    std::memcpy(&value, data, sizeof(value));
                    v = value;
                }
                else if (mat.type() == CV_64FC1){
                    double value = 0.0f;
                    std::memcpy(&value, data, sizeof(value));
                    v = value;
                }
                else if (mat.type() == CV_32SC1){
                    int32_t value = 0;
                    std::memcpy(&value, data, sizeof(value));
                    v = value;
                }
                else if (mat.type() == CV_64FC1){
                    int64_t value = 0;
                    std::memcpy(&value, data, sizeof(value));
                    v = value;
                }
                else if (mat.type() == CV_8UC1){
                    uint8_t value = 0;
                    std::memcpy(&value, data, sizeof(value));
                    v = value;
                }
                else{
                    assert(false && "element type is not supported here!");
                }
                if (v != 0.0)
                    triplets.emplace_back(jj, ii, v); // col - row - val
            }

            // make triplets ordered in column indices to satisfy matlab interface
            std::sort(triplets.begin(), triplets.end());

            // fill in matlab data
            auto sr = mxGetPr(ma);
            //auto si = mxGetPi(ma);
            auto irs = mxGetIr(ma);
            auto jcs = mxGetJc(ma);
            jcs = (mwIndex*)mxRealloc(jcs, (im.size(1) + 1) * sizeof(mwIndex));
            std::fill(jcs, jcs + (im.size(1) + 1), 0);
            mxSetJc(ma, jcs);

            for (int i = 0; i < triplets.size(); i++){
                sr[i] = std::get<2>(triplets[i]);
                int ii = std::get<1>(triplets[i]) + 1;
                int jj = std::get<0>(triplets[i]) + 1;
                irs[i] = ii - 1;
                jcs[jj] ++;
            }

            for (int j = 1; j < (im.size(1) + 1); j++){
                jcs[j] += jcs[j - 1];
            }

            return MXArray(ma, dos);
        }





        MAT::MAT() : _fp(nullptr) {}
        MAT::MAT(const std::string & fname, const std::string & mode) : _fname(fname){
            _fp = matOpen(fname.c_str(), mode.c_str());
        }
        MAT::MAT(const std::string & fname, OpeningMode mode) : _fname(fname), _fp(nullptr) {
            std::string m;
            switch (mode){
            case Read:
                m = "r";
                break;
            case Update:
                m = "u";
                break;
            case Write:
                m = "w";
                break;
            case Write_4:
                m = "w4";
                break;
            case Write_CompressedData:
                m = "wz";
                break;
            case Write_7_3:
                m = "w7.3";
                break;
            default:
                return;
            }
            _fp = matOpen(fname.c_str(), m.c_str());
        }

        MAT::MAT(MAT && a) {
            _fname = a._fname;
            a._fname.clear();
            _fp = a._fp;
            a._fp = nullptr;
        }
        MAT & MAT::operator = (MAT && a) {
            std::swap(_fname, a._fname);
            std::swap(_fp, a._fp);
            return *this;
        }


        MAT::~MAT() {
            if (_fp){
                matClose(_fp);
                _fp = nullptr;
            }
        }


        std::vector<std::string> MAT::varNames() const{
            int num = 0;
            char ** names = matGetDir(_fp, &num);
            if (!names)
                return std::vector<std::string>();
            std::vector<std::string> vnames(num);
            for (int i = 0; i < num; i++){
                vnames[i] = names[i];
            }
            mxFree(names);
            return vnames;
        }


    }

}
