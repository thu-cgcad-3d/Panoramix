#ifdef USE_MATLAB
#include <mex.h>
#include <engine.h>
#include "matlab_api.hpp"
#endif

#include "matlab_engine.hpp"

namespace panoramix {
    namespace misc {

        using namespace core;

#ifdef USE_MATLAB

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

        class Engine {
        public:
            inline explicit Engine(const char * cmd = nullptr, int bufferSize = 5e4)
                : _engine(nullptr), _buffer(nullptr) {
                _engine = engOpen(cmd);
                if (_engine){
                    engSetVisible(_engine, false);
                    _buffer = new char[bufferSize];
                    std::memset(_buffer, 0, bufferSize);
                    engOutputBuffer(_engine, _buffer, bufferSize);
                    std::cout << "MatlabEngine Engine Launched" << std::endl;
                    engEvalString(_engine, ("cd " + MatlabEngine::DefaultCodeDir() + "; startup; pwd").c_str());
                    std::cout << _buffer << std::endl;
                }
            }
            inline ~Engine() {
                if (_engine){
                    engClose(_engine);
                    std::cout << "MatlabEngine Engine Closed" << std::endl;
                }
                delete[] _buffer;
            }
            inline engine * eng() const { return _engine; }
            inline char * buffer() const { return _buffer; }
        private:
            engine * _engine;
            char * _buffer;
        };

        static Engine * _Engine = nullptr;
        static int _EngineRefCount = 0;


        bool MatlabEngine::IsBuilt() { return true; }

        bool MatlabEngine::Start(){
            _EngineRefCount++;
            if (_EngineRefCount > 0 && !_Engine){
                _Engine = new Engine;
                if (_Engine->eng() == nullptr){
                    delete _Engine;
                    _Engine = nullptr;
                }
            }
            std::cout << "MatlabEngine Engine RefCount: " << _EngineRefCount << std::endl;
            return _Engine != nullptr;
        }

        void MatlabEngine::Close(){
            if (_EngineRefCount <= 0){
                _EngineRefCount = 0;
                return;
            }
            _EngineRefCount--;
            std::cout << "MatlabEngine Engine RefCount: " << _EngineRefCount << std::endl;
            if (_EngineRefCount <= 0){
                _EngineRefCount = 0;
                if (_Engine){
                    delete _Engine;
                    _Engine = nullptr;
                }
            }
        }

        bool MatlabEngine::Started(){
            return _Engine && _Engine->eng();
        }

        

        std::string MatlabEngine::DefaultCodeDir() { return MATLAB_CODE_DIR; }


        bool MatlabEngine::RunScript(const char * cmd) {
            if (!_Engine)
                return false;
            bool ret = engEvalString(_Engine->eng(), cmd) == 0;
            if (strlen(_Engine->buffer()) > 0){
                std::cout << "[Message when executing '" << cmd << "']:\n" << _Engine->buffer() << std::endl;
            }
            return ret;
        }


        const char * MatlabEngine::LastMessage(){
            return _Engine->buffer();
        }

        bool MatlabEngine::PutVariable(const char * name, CVInputArray mat){
            if (!_Engine)
                return false;
            
            MXA ma(mat, true);
            if (!ma)
                return false;

            int result = engPutVariable(_Engine->eng(), name, ma);
            return result == 0;
        }


        bool MatlabEngine::GetVariable(const char * name, CVOutputArray mat, bool lastDimIsChannel){
            if (!_Engine)
                return false;

            if (!mat.needed())
                return true;
            
            MXA ma(engGetVariable(_Engine->eng(), name), true);
            if (ma.null())
                return false;

            bool result = ma.toCVOutputArray(mat, lastDimIsChannel);
            return result;
        }

        bool MatlabEngine::GetVariable(const char * name, double & a){
            if (!_Engine)
                return false;

            MXA ma(engGetVariable(_Engine->eng(), name), true);
            if (ma.null())
                return false;

            a = ma;
            return true;
        }


        bool MatlabEngine::GetVariable(const char * name, std::string & a){
            if (!_Engine)
                return false;

            MXA ma(engGetVariable(_Engine->eng(), name), true);
            if (ma.null())
                return false;

            a = ma;
            return true;
        }


        bool MatlabEngine::PutVariable(const char * name, const cv::SparseMat & mat){
            if (!_Engine)
                return false;
            
            MXA ma(mat, true);
            if (ma.null())
                return false;

            int result = engPutVariable(_Engine->eng(), name, ma);
            return result == 0;
        }


#else
        bool MatlabEngine::IsBuilt() {return false;}
        bool MatlabEngine::Start() {return false;}
        void MatlabEngine::Close() {}
        bool MatlabEngine::Started(){ return false;}
        std::string MatlabEngine::DefaultCodeDir() { return ""; }
        bool MatlabEngine::RunScript(const char *) {return false;}
        const char * MatlabEngine::LastMessage() {return nullptr;}
        bool MatlabEngine::PutVariable(const char * name, CVInputArray a) {return false;}
        bool MatlabEngine::GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel) {return false;}
        bool MatlabEngine::PutVariable(const char * name, const cv::SparseMat & mat) { return false; }
#endif

    }
}