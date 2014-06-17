#include "regions_net.hpp"

namespace panoramix {
    namespace rec {

        RegionsNet::RegionsNet(const Image & image, const Params & params) 
            : _image(image), _params(params)/*, _regionsRTree(RegionDataBoundingBoxFunctor(_regions))*/{}

        namespace {

            void ComputeRegionProperties(const Image & regionMask, const Image & image,
                Vec2 & center, double & area, Box2 & boundingBox) {
                assert(regionMask.depth() == CV_8U && regionMask.channels() == 1);
                center = { 0, 0 };
                area = 0;
                boundingBox = Box2();
                for (auto i = regionMask.begin<uint8_t>(); i != regionMask.end<uint8_t>(); ++i){
                    if (*i) { // masked
                        area += 1.0;
                        PixelLoc p = i.pos();
                        boundingBox |= BoundingBox(p);
                        center += Vec2(p.x, p.y);
                    }
                }
                center /= area;
            }

        }

        void RegionsNet::buildNetAndComputeGeometricFeatures() {
            _segmentedRegions = _params.segmenter(_image);
            int regionNum = static_cast<int>(MinMaxValOfImage(_segmentedRegions).second) + 1;
            _regions.internalElements<0>().reserve(regionNum);
            for (int i = 0; i < regionNum; i++){
                RegionData vd;
                vd.borderLength = 0.0;
                vd.regionMask = (_segmentedRegions == i);
                ComputeRegionProperties(vd.regionMask, _image, vd.center, vd.area, vd.boundingBox);

                // find contour of the region
                cv::Mat regionMaskCopy;
                vd.regionMask.copyTo(regionMaskCopy);
                std::vector<std::vector<PixelLoc>> contours;
                cv::findContours(regionMaskCopy, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
                vd.contour = contours.front();
                assert(!vd.contour.empty() && "no contour? impossible~");
               
                auto vh = _regions.add(vd);
                //_regionsRTree.insert(vh);  
            }

            // find connections
            int width = _segmentedRegions.cols;
            int height = _segmentedRegions.rows;
            std::vector<double> boundaryLengths(regionNum * regionNum, 0.0);
            std::set<std::pair<int, int>> connections;
            for (int y = 0; y < height-1; y++){
                for (int x = 0; x < width-1; x++){
                    PixelLoc p1(x, y), p2(x + 1, y), p3(x, y + 1), p4(x + 1, y + 1);
                    int r1 = _segmentedRegions.at<int32_t>(p1),
                        r2 = _segmentedRegions.at<int32_t>(p2),
                        r3 = _segmentedRegions.at<int32_t>(p3),
                        r4 = _segmentedRegions.at<int32_t>(p4);
                    if (r1 != r2){
                        double blen = 1.0;
                        boundaryLengths[r1 * regionNum + r2] += blen;
                        boundaryLengths[r2 * regionNum + r1] += blen;
                        _regions.data(RegionHandle(r1)).borderLength += blen;
                        _regions.data(RegionHandle(r2)).borderLength += blen;
                        connections.insert(std::make_pair(r1, r2));
                        connections.insert(std::make_pair(r2, r1));
                    }
                    if (r1 != r3){
                        double blen = 1.0;
                        boundaryLengths[r1 * regionNum + r3] += blen;
                        boundaryLengths[r3 * regionNum + r1] += blen;
                        _regions.data(RegionHandle(r1)).borderLength += blen;
                        _regions.data(RegionHandle(r3)).borderLength += blen;
                        connections.insert(std::make_pair(r1, r3));
                        connections.insert(std::make_pair(r3, r1));
                    }
                    if (r1 != r4){
                        double blen = 1.4;
                        boundaryLengths[r1 * regionNum + r4] += blen;
                        boundaryLengths[r4 * regionNum + r1] += blen;
                        _regions.data(RegionHandle(r1)).borderLength += blen;
                        _regions.data(RegionHandle(r4)).borderLength += blen;
                        connections.insert(std::make_pair(r1, r4));
                        connections.insert(std::make_pair(r4, r1));
                    }
                }
            }

            for (auto con : connections) {
                if (con.first < con.second){
                    BoundaryData hd;
                    hd.length = boundaryLengths[con.first * regionNum + con.second];
                    _regions.add<1>({ RegionHandle(con.first), RegionHandle(con.second) }, hd);
                }
            }

            // todo: compute boundary line, boundary roughness
            //NOT_IMPLEMENTED_YET();
        }

        void RegionsNet::computeImageFeatures() {
            // NOT_IMPLEMENTED_YET();
        }

    }
}