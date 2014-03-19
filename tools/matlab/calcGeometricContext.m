function calcGeometricContext(im)

eclassifier = {};
labelclassifier = {};
segclassifier = {};
ecal = [];

load Classifiers_gc;
	
tmpimdir = [tempdir 'temp_name.ppm'];
tmpsegdir = [tempdir 'temp_seg.pnm'];
imwrite(im, tmpimdir);
system(['"segmentation" 0.8 100 100 ' tmpimdir ' ' tmpsegdir]);
imseg = processSuperpixelImage(tmpsegdir);
delete(tmpimdir);
delete(tmpsegdir);

seg_im = label2rgb(imseg.segimage, 'jet', 'w', 'shuffle');

% initial gc
spfea = mcmcGetSuperpixelData(im2double(im), imseg);
[edgefea, adjlist, ~, ~] = mcmcGetEdgeData(imseg, spfea);	
confidences = test_boosted_dt_mc(eclassifier, edgefea);
confidences = 1 ./ (1+exp(ecal(1)*confidences+ecal(2)));

smaps = msCreateMultipleSegmentations(confidences, adjlist, ...
	imseg.nseg, [5 15 25 35 40 60 80 100]);
imdata = mcmcComputeImageData(im2double(im), imseg);
   segfea = {};
for jdx=1:size(smaps, 2)
	if max(smaps(:, jdx)) > 0
		segfea{jdx} = mcmcGetSegmentFeatures(imseg, ...
			spfea, imdata, smaps(:, jdx), 1:max(smaps(:, jdx)));
	end
end

initSlabelConfidence = msTest(imseg, segfea, {smaps}, labelclassifier, segclassifier, 1);
initSlabelConfidence = initSlabelConfidence{1};
[slabelConfMap, ~, slabelConfMap_im] = ...
	composeSegLabel(imseg.segimage, initSlabelConfidence);
slabelConfMap_im = slabelConfMap_im / max(slabelConfMap_im(:));
imshow(slabelConfMap_im);

end