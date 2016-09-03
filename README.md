# Panoramix
Codes of the CVPR'16 paper: 
<cite>
Hao Yang and Hui Zhang. "Efficient 3D Room Shape Recovery From a Single Panorama." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2016.
</cite>

Author: Hao Yang (yangh2007@gmail.com)

## Build

### Requirements
* Only Visual Studio 2015 is supported currently.
* OpenCV is installed;
* Qt5 is installed;
* MATLAB is installed;
* CVX is installed;
* The [MATLABTools](https://github.com/YANG-H/MATLABTools) repo is downloaded, the CVX path is set in file `startup.m`;
* CMake is installed.

### Build with CMake
* Set `OpenCV_DIR` to the path of OpenCV;
* Set `Qt_DIR` to the path of Qt;
* Set `MATLAB_CODE_DIR` to the path of MATLABTools;
* Build.