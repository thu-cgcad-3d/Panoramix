#ifndef PANORAMIX_REC_VIEWS_NET_HPP
#define PANORAMIX_REC_VIEWS_NET_HPP

#include "lines_net.hpp"
#include "regions_net.hpp"
 
namespace panoramix {
    namespace rec {

        using namespace core;

        // net of line segments
        class ViewsNet {
        public:
            struct Params {
                Params();
                
            };

            struct ViewData {
                
            };
            struct ViewConnectionData {
                
            };
            using ViewsGraph = GraphicalModel02 < ViewData, ViewConnectionData >;
            using ViewHandle = HandleAtLevel < 0 > ;
            using ViewConnectionHandle = HandleAtLevel < 1 > ;
        
        public:
            explicit ViewsNet(const Params & params = Params());
            inline const ViewsGraph & views() const {return _views;}
            
        private:
            ViewsGraph _views;
            Params _params;
        };

    }
}
 
#endif