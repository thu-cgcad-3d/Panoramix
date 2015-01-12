#ifdef USE_MATLAB
#include <engine.h>
#endif

#include "matlab.hpp"

namespace panoramix {
    namespace core {

#ifdef USE_MATLAB
        namespace {

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
        }

        class MatlabEngine {
        public:
            inline explicit MatlabEngine(const char * cmd = nullptr, int bufferSize = 1024) 
                : _engine(nullptr), _buffer(nullptr) { 
                _engine = engOpen(cmd);
                if (_engine){
                    engSetVisible(_engine, false);
                    _buffer = new char[bufferSize];
                    std::memset(_buffer, 0, bufferSize);
                    engOutputBuffer(_engine, _buffer, bufferSize);
                    std::cout << "Matlab Engine Launched" << std::endl;
                    engEvalString(_engine, ("cd " + Matlab::DefaultCodeDir() + "; startup; pwd").c_str());
                    std::cout << _buffer << std::endl;
                }
            }
            inline ~MatlabEngine() { 
                if (_engine){
                    engClose(_engine);
                    std::cout << "Matlab Engine Closed" << std::endl;
                }
                delete[] _buffer;
            }
            inline engine * eng() const { return _engine; }
            inline char * buffer() const { return _buffer; }
        private:
            engine * _engine;
            char * _buffer;
        };
        
        static MatlabEngine _Engine;

        bool Matlab::IsBuilt() { return true; }

        std::string Matlab::DefaultCodeDir() { return MATLAB_CODE_DIR; }

        bool Matlab::IsUsable() {
            return _Engine.eng() != nullptr;
        }

        bool Matlab::RunScript(const char * cmd) {
            if (!_Engine.eng())
                return false;
            bool ret = engEvalString(_Engine.eng(), cmd) == 0;
            if (strlen(_Engine.buffer()) > 0){
                std::cout << "[Message when executing '" << cmd << "']:\n" << _Engine.buffer() << std::endl;
            }
            return ret;
        }

        const char * Matlab::LastMessage(){
            return _Engine.buffer();
        }

        bool Matlab::PutVariable(const char * name, CVInputArray mat){
            if (!_Engine.eng())
                return false;
            
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
                return false;

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

            int result = engPutVariable(_Engine.eng(), name, ma);
            mxDestroyArray(ma);
            return result == 0;
        }

        bool Matlab::GetVariable(const char * name, CVOutputArray mat, bool lastDimIsChannel){
            if (!_Engine.eng())
                return false;

            if (!mat.needed())
                return true;
            
            mxArray * ma = engGetVariable(_Engine.eng(), name);
            if (!ma)
                return false;

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
                mxDestroyArray(ma);
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

            mxDestroyArray(ma);
            return true;
        }


#else
        bool Matlab::IsBuilt() { return false; }
        std::string Matlab::DefaultCodeDir() { return ""; }
        bool Matlab::IsUsable() {return false;}
        bool Matlab::RunScript(const char *) {return false;}
        const char * Matlab::LastMessage() {return nullptr;}
        bool Matlab::PutVariable(const char * name, CVInputArray a) {return false;}
        bool Matlab::GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel) {return false;}
#endif

    }
}