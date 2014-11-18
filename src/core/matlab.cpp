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
            inline explicit MatlabEngine(const char * cmd = nullptr) { 
                _engine = engOpen(cmd);
                engSetVisible(_engine, false);
            }
            inline ~MatlabEngine() { 
                engClose(_engine); 
            }
            inline engine * eng() const { return _engine; }
        private:
            engine * _engine;
        };
        
        static MatlabEngine _Engine;

        bool Matlab::IsBuilt() { return true; }

        bool Matlab::IsUsable() {
            return _Engine.eng() != nullptr;
        }

        bool Matlab::RunScript(const char * cmd) {
            return engEvalString(_Engine.eng(), cmd) == 0;
        }

        bool Matlab::PutVariable(const char * name, const Image & im){
            int channelNum = im.channels();
            // collect all dimensions of im
            std::vector<mwSize> dims(im.dims + 1);
            for (int i = 0; i < im.dims; i++)
                dims[i] = im.size[i];
            dims.back() = channelNum;
            // create mxArray
            mxArray* ma = mxCreateNumericArray(dims.size(), dims.data(), CVDepthToMxClassID(im.depth()), mxREAL);
            if (!ma)
                return false;
            uint8_t * mad = (uint8_t*)mxGetData(ma);
            const size_t szForEachElem = im.elemSize1();

            //for (auto i = im.begin())

            for (mwIndex i = 0; i < im.rows; i++){
                for (mwIndex j = 0; j < im.cols; j++){
                    for (mwIndex k = 0; k < channelNum; k++) {
                        mwIndex idx[] = { i, j, k };
                        const uint8_t * fromDataHead = im.ptr(i, j) + k * szForEachElem;
                        uint8_t* toDataHead = mad + mxCalcSingleSubscript(ma, 3, idx) * szForEachElem;
                        std::memcpy(toDataHead, fromDataHead, szForEachElem);
                    }
                }
            }
            int result = engPutVariable(_Engine.eng(), name, ma);
            mxDestroyArray(ma);
            return result == 0;
        }

        bool Matlab::GetVariable(const char * name, Image & im){
            mxArray * ma = engGetVariable(_Engine.eng(), name);
            if (!ma)
                return false;
            int d = mxGetNumberOfDimensions(ma);
            assert((d == 2 || d == 3) && "dimension num of the variable must be 2 or 3 for representing as an Image");
            const mwSize * dims = mxGetDimensions(ma);
            int rows = dims[0];
            int cols = dims[1];
            int channels = d == 3 ? dims[2] : 1;
            const size_t szForEachElem = mxGetElementSize(ma);
            
            int depth = MxClassIDToCVDepth(mxGetClassID(ma));
            if (depth == -1){
                mxDestroyArray(ma);
                return false;
            }

            const uint8_t * mad = (const uint8_t*)mxGetData(ma);
            im.create(rows, cols, CV_MAKETYPE(depth, channels));
            for (mwIndex i = 0; i < rows; i++){
                for (mwIndex j = 0; j < cols; j++){
                    for (mwIndex k = 0; k < channels; k++) {
                        mwIndex idx[] = { i, j, k };
                        uint8_t * toDataHead = im.ptr(i, j) + k * szForEachElem;
                        const uint8_t* fromDataHead = mad + mxCalcSingleSubscript(ma, 3, idx) * szForEachElem;
                        std::memcpy(toDataHead, fromDataHead, szForEachElem);
                    }
                }
            }
            mxDestroyArray(ma);
            return true;
        }


#else
        bool Matlab::IsBuilt() { return false; }
        bool Matlab::IsUsable() {return false;}
        bool Matlab::RunScript(const char *) {return false;}
        bool Matlab::PutVariable(const char * name, const Image & im) {return false;}
        bool Matlab::GetVariable(const char * name, Image & im) {return false;}
#endif

    }
}