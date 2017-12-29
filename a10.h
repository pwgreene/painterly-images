#include "matrix.h"
#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include "assert.h"
#include "image_io.h"
#include <Halide.h>
#include <stdlib.h>
#include <time.h>


#ifndef A10_H_PHUDVTKB
#define A10_H_PHUDVTKB

// Write your declarations here, or extend the Makefile if you add source
// files

//im: the image to draw to. y,x: where to draw in out. color: the color of the stroke. texture: the texture of the stroke.
void brush(Image &output, int x, int y, Vec3f color, const Image &texture);

//Paints with all brushed at the same scale using importance sampling.
void singleScalePaint(const Image &im, Image &output, const Image &importance, const Image &texture, int N=1000, int size=10, float noise=0.3);

//First paints at a coarse scale using all 1's for importance sampling, then paints again at size/4 scale using the sharpness map for importance sampling.
Image painterly(Image &im, const Image &texture, int N=10000, int size=50, float noise=0.3);
Halide::Func painterlyHalide(const Halide::Tools::Image<float> &im, const Halide::Tools::Image<float> &texture, int N=10000, int size=50, float noise=0.3);

//Return an image that holds the angle of the smallest eigenvector of the structure tensor at each pixel. If you have a 3 channel image as input, just set all three channels to be the same value theta.
Image computeAngles(const Image &im);
Halide::Func computeAnglesHalide(Halide::Func blurLumi, Halide::Tools::Image<float> gaussKernel, int w, int h, int radius);

//same as single scale paint but now the brush strokes will be oriented according to the angles in thetas.
void singleScaleOrientedPaint(const Image &im, Image &output, const Image &thetas, const Image &importance, const Image &texture, int N=10000, int size=10, float noise=0.3, int nAngles=36);

//same as painterly but computes and uses the local orientation information to orient strokes.
Image orientedPaint(const Image &im, const Image &texture, int N=7000, int size=50, float noise=0.3);

Image sharpnessMap(const Image &im, float sigma=1);
Image computeTensor(const Image &im, float sigmaG=1, float factorSigma=4);
//draw white brushes rotated to each angle (for debugging purposes)
Image drawAngles(const Image &angles, const Image &texture, int size=10);
vector<Image> rotateBrushes(const Image &texture, int nAngles);
Halide::Func rotateBrushesHalide (Halide::Func texture, int nAngles=36);



#endif /* end of include guard: A10_H_PHUDVTKB */

