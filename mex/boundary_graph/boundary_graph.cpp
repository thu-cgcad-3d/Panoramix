#include "mex.h"

#include "../../src/core/utilities.hpp"

#include "boundary_graph.hpp"

namespace mex {

    BoundaryGraph::BoundaryGraph(const Image & im, const Imaged7 & gc, const Imaged & d) {

        image = im;
        geometricContext = gc;
        depth = d;

        int segsNum = 0;
        Imagei segs;
        SegmentationExtractor segmenter;
        segmenter.params().algorithm = SegmentationExtractor::GraphCut;
        std::tie(segs, segsNum) = segmenter(image);

        std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
        FindContoursOfRegionsAndBoundaries(segs, segsNum, boundaryEdges, 3);

        graph.internalElements<0>().reserve(segsNum);
        graph.internalElements<1>().reserve(boundaryEdges.size());

        //#pragma parallel for
        for (int i = 0; i < segsNum; i++){
            Image regionMask = (segs == i);
            Region rd;
            // find contour of the region
            // CV_RETR_EXTERNAL: get only the outer contours
            cv::findContours(regionMask, rd.contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
            if (rd.contours.empty()){
                continue;
            }
            rd.perimeter = 0.0;
            for (auto & cs : rd.contours){
                for (int k = 1; k < cs.size(); k++){
                    rd.perimeter += Distance(cs[k - 1], cs[k]);
                }
            }
            rd.meanColor = Vec3(0, 0, 0);
            rd.area = 0;

            // colorhist TODO

            // plane TODO
            if (!d.empty()){

            }
            graph.add(std::move(rd));
        }

        for (auto it = segs.begin(); it != segs.end(); ++it){
            int segid = *it;
            auto & rd = graph.internalElements<0>()[segid];

            rd.data.area++;
            rd.data.center += Point2(it.pos().x, it.pos().y);
            rd.data.meanColor += (Vec3)gui::ColorFromImage(im, it.pos());
            
            if (!gc.empty()){
                rd.data.meanGC += gc(it.pos());
            }
        }

        for (int i = 0; i < segsNum; i++){
            auto & rd = graph.internalElements<0>()[i];
            rd.data.center /= rd.data.area;
            rd.data.meanColor /= rd.data.area;
            rd.data.meanGC /= rd.data.area;
        }

        for (auto & bep : boundaryEdges){
            auto & rids = bep.first;
            Boundary bd;
            // todo
            bd.pixels = std::move(bep.second);
            graph.add<1>({ RegionHandle(rids.first), RegionHandle(rids.second) }, std::move(bd));
        }

        segmentedImage = segs;
        segmentsNum = segsNum;

    }

    void BoundaryGraph::print() const{
        printf("BoundaryGraph with %d regions and %d boundaries\n", 
            graph.internalElements<0>().size(), graph.internalElements<1>().size());
    }

    Imaged BoundaryGraph::computeFeaturesMat() const{
        NOT_IMPLEMENTED_YET();
    }

    Imaged BoundaryGraph::computeLabelsVec() const{
        NOT_IMPLEMENTED_YET();
    }

    void BoundaryGraph::installLabelsVec(const Imaged & labels){
        NOT_IMPLEMENTED_YET();
    }

}