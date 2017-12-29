#include <iostream>
#include "a10.h"

using namespace std;
using namespace Halide;

// Write your implementations here, or extend the Makefile if you add source
// files

void brush(Image &output, int x, int y, Vec3f color, const Image &texture) {
    //make sure we have a square texture
    int halfTexWidth = texture.width()/2;
    int halfTexHeight = texture.height()/2;
    if (x < halfTexWidth  || x > output.width()-halfTexWidth ||
        y < halfTexHeight || y > output.height()-halfTexHeight) {
        return; //do nothing if outside valid area
    }
    
    float alpha;
    for (int i = -halfTexWidth; i < halfTexWidth; i++) {
        for (int j = -halfTexHeight; j < halfTexHeight; j++) {
            for (int c = 0; c < output.channels(); c++) {
                alpha = texture(i+halfTexWidth, j+halfTexHeight);
//                cout << alpha << "\n";
                output(x+i, y+j, c) = output(x+i, y+j, c)*(1-alpha) + color[c]*alpha;
            }
        }
    }
}

void singleScalePaint(const Image &im, Image &output, const Image &importance, const Image &texture, int N, int size, float noise) {
    
    //rescale texture to maximum size
    int maxDim = max(texture.width(), texture.height());
    Image scaledTexture = scaleNN(texture, (float)(size)/maxDim);
    
    float avgImportance = importance.mean();
    //rescale N
    N /= avgImportance;
    
    int randX, randY;
    Vec3f color;
    float modulation;
    for (int i = 0; i < N; i++) {
        randX = rand() % output.width();
        randY = rand() % output.height();
        //rejection sampling
        if ((float)(rand())/RAND_MAX < importance(randX, randY)) {
            modulation = 1-noise/2 + noise*((float)(rand())/RAND_MAX);
            for (int c = 0; c < output.channels(); c++) {
                color[c] = im(randX, randY, c) * modulation;
            }
            brush(output, randX, randY, color, scaledTexture);
        }
    }
}

Image painterly(Image &im, const Image &texture, int N, int size, float noise) {
    Image output(im.width(), im.height(), im.channels());
    output.set_color(0, 0, 0);
    
    Image importanceConst(im.width(), im.height());
    importanceConst.set_color(1);
    Image importanceSharpness = sharpnessMap(im);
    
    //first pass
    singleScalePaint(im, output, importanceConst, texture, N, size, noise);
    //second pass, finer detail
    singleScalePaint(im, output, importanceSharpness, texture, N, size/4, noise);
    
    return output;
}

Image computeAngles(const Image &im) {
    Image output(im.width(), im.height(), 1);
    Image tensor = computeTensor(im);
//    Matrix tensor(2, 2);
//    Eigen::EigenSolver<Matrix> es;
    float trace, det, eval1, eval2, vecx, vecy;
    for (int x = 0; x < tensor.width(); x++) {
        for (int y = 0; y < tensor.height(); y++) {
            trace = tensor(x, y, 0) + tensor(x, y, 2);
            det = tensor(x, y, 0)*tensor(x, y, 2) - tensor(x, y, 1)*tensor(x, y, 1);
            eval1 = trace*.5f + pow((pow(trace, 2)*.25f - det), .5f);
            eval2 = trace*.5f - pow((pow(trace, 2)*.25f - det), .5f);
            vecx = eval1 < eval2 ? eval1-tensor(x, y, 2) : eval2-tensor(x, y, 2);
            vecy = tensor(x, y, 1);
            output(x, y, 0) = atan2(vecy, vecx) + 3.1415926f;
//            tensor(0, 0) = tensorIm(x, y, 0); tensor(0, 1) = tensorIm(x, y, 1);
//            tensor(1, 0) = tensorIm(x, y, 1); tensor(1, 1) = tensorIm(x, y, 2);
//            es = Eigen::EigenSolver<Matrix>(tensor);
//            //get smallest eigenvector
//             if (abs(es.eigenvalues()[0].real()) < abs(es.eigenvalues()[1].real()))
//            {
//                output(x, y, 0) = atan2(es.eigenvectors().col(0)[1].real(),
//                                        es.eigenvectors().col(0)[0].real()) + M_PI;
//            } else {
//                output(x, y, 0) = atan2(es.eigenvectors().col(1)[1].real(),
//                                        es.eigenvectors().col(1)[0].real()) + M_PI;
//            }
        }
    }
    return output;
}

void singleScaleOrientedPaint(const Image &im, Image &output, const Image &thetas, const Image &importance, const Image &texture, int N, int size, float noise, int nAngles) {

    //rescale texture to maximum size
    int maxDim = max(texture.width(), texture.height());
    Image scaledTexture = scaleNN(texture, (float)(size)/maxDim);
    Image rotatedTexture(scaledTexture.width(), scaledTexture.height(), scaledTexture.channels());
    vector<Image> brushes = rotateBrushes(scaledTexture, nAngles);
    
    float avgImportance = importance.mean();
    //rescale N
    N /= avgImportance;
    
    int randX, randY;
    Vec3f color;
    float modulation;
    for (int i = 0; i < N; i++) {
        randX = rand() % output.width();
        randY = rand() % output.height();
        //rejection sampling
        if ((float)(rand())/RAND_MAX < importance(randX, randY)) {
//            cout << (int)(thetas(randX, randY, 0)*nAngles/(2*M_PI)) % nAngles << "\n";
            rotatedTexture = brushes[(int)(thetas(randX, randY, 0)*nAngles/(2*M_PI)) % nAngles];
            modulation = 1-noise/2 + noise*((float)(rand())/RAND_MAX);
            for (int c = 0; c < output.channels(); c++) {
                color[c] = im(randX, randY, c) * modulation;
            }
            brush(output, randX, randY, color, rotatedTexture);
        }
    }
}

Image orientedPaint(const Image &im, const Image &texture, int N, int size, float noise) {
    Image output(im.width(), im.height(), im.channels());
    output.set_color(0, 0, 0);
    
    Image importanceConst(im.width(), im.height());
    importanceConst.set_color(1);
    Image importanceSharpness = sharpnessMap(im);
    Image thetas = computeAngles(im);
    
    //first pass
    singleScaleOrientedPaint(im, output, thetas, importanceConst, texture, N, size, noise);
    //second pass, finer detail
    singleScaleOrientedPaint(im, output, thetas, importanceSharpness, texture, N, size/4, noise);
    
    return output;
}

//helper functions

Image sharpnessMap(const Image &im, float sigma) {
    Image lumi = lumiChromi(im)[0];
    Image blur = gaussianBlur_separable(lumi, sigma);
    Image high = lumi-blur;
    Image energy = high*high;
    Image sharpness=gaussianBlur_separable(energy, 4*sigma);
    sharpness = sharpness/sharpness.max();
    return sharpness;
}

Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
    // Compute xx/xy/yy Tensor of an image. (stored in that order)
    Image output(im.width(), im.height(), 3);
    Image lumi_blurred = gaussianBlur_separable(lumiChromi(im)[0], sigmaG);
    Image gradx = gradientX(im);
    Image grady = gradientY(im);
    for (int x = 0; x < output.width(); x++) {
        for (int y = 0; y < output.height(); y++) {
            output(x, y, 0) = gradx(x, y, 0) * gradx(x, y, 0);
            output(x, y, 1) = gradx(x, y, 0) * grady(x, y, 0);
            output(x, y, 2) = grady(x, y, 0) * grady(x, y, 0);
        }
    }
    return gaussianBlur_separable(output, sigmaG*factorSigma);
}

Image drawAngles(const Image &angles, const Image &texture, int size) {
    Image output(angles.width(), angles.height(), angles.channels());
    int maxDim = max(texture.width(), texture.height());
    Image scaledTexture = scaleLin(texture, (float)(size)/maxDim);
    output.set_color(0, 0, 0);
    Vec3f white(1, 1, 1);
    Image rotatedTexture(texture.width(), texture.height(), 1);
    for (int x = 0; x < output.width(); x+=size) {
        for (int y = 0; y < output.height(); y+=size) {
//            cout << x << ", " << y << "\n";
            rotatedTexture = rotate(scaledTexture, angles(x, y, 0));
            brush(output, x, y, white, rotatedTexture);
        }
    }
    return output;
}

vector<Image> rotateBrushes(const Image &texture, int nAngles) {
    vector<Image> brushes;
    float theta;
    for (int i = 0; i < nAngles; i++) {
        theta = 2*M_PI*i/nAngles;
        brushes.push_back(rotate(texture, -theta));
    }
    return brushes;
}



//Halide versions

Func rotateBrushesHalide (Func texture, int nAngles, int w, int h) {
    Func brushes("brushes");
    Func texClamped("texClamped");
    Var x("x"), y("y"), i("i");
    Var xo("xo"), yo("yo"), xi("xi"), yi("yi"), io("io"), ii("ii");
    int centerX = (w-1.f)/2.f;
    int centerY = (h-1.f)/2.f;
    texClamped(x, y) = texture(clamp(x, 0, w-1), clamp(y, 0, h-1));

    //rotate the texture
    Expr theta = -2*3.1415926f*i / nAngles;
    Expr xR = (x - centerX)*cos(theta) + (centerY - y)*sin(theta) + centerX;
    Expr yR = centerY - (-(x - centerX)*sin(theta) + (centerY - y)*cos(theta));
    //linear interpolation, see basicImageManipulation.cpp for reference c++ implementation
    Expr xf = cast<int>(floor(xR));
    Expr yf = cast<int>(floor(yR));
    Expr xc = xf+1;
    Expr yc = yf+1;
    Expr yalpha = yR - yf;
    Expr xalpha = xR - xf;
    Expr tl = texClamped(xf, yf); Expr tr = texClamped(xc, yf);
    Expr bl = texClamped(xf, yc); Expr br = texClamped(xc, yc);
    Expr topL = tr*xalpha + tl*(1.f - xalpha);
    Expr botL = br*xalpha + bl*(1.f - xalpha);
    Expr retv = botL*yalpha + topL*(1.f - yalpha);
    brushes(x, y, i) = retv;
    brushes.parallel(i);
    brushes.compute_root();
    return brushes;
    }

    void apply_compute_root(Func F) {
    map<string,Internal::Function> flist = Internal::find_transitive_calls(F.function());
    flist.insert(std::make_pair(F.name(), F.function()));
    map<string,Internal::Function>::iterator fit;
    for (fit=flist.begin(); fit!=flist.end(); fit++) {
        Func f(fit->second);
        f.compute_root();
        // cout << "Warning: applying default schedule to " << f.name() << endl;
    }
    cout << endl;
    }

Func computeAnglesHalide(Halide::Func blurLumi, Halide::Tools::Image<float> gaussKernel, int w, int h, int radius) {

    int fwidth = 2 * radius + 1;
    RDom rx(0, fwidth);
    Func tensorClamped;
    Func blurTensorX;
    Func clampedBlurTensorX;
    Func blurTensor;
    Func output;

    //compute gradX and gradY
    Func gradX, gradY, clampedBlurLumi;
    Var x, y, c, xo, yo, xi, yi;
    Func tensor;
    Expr x_clamped = clamp(x, 0, w-1);
    Expr y_clamped = clamp(y, 0, h-1);
    clampedBlurLumi(x, y) = cast<float>(blurLumi(x_clamped, y_clamped));
    gradX(x, y) =cast<float> (  -clampedBlurLumi(x-1, y-1) + clampedBlurLumi(x+1, y-1)
                              + -clampedBlurLumi(x-1, y)*2 + 2*clampedBlurLumi(x+1, y)
                              + -clampedBlurLumi(x-1, y+1) + clampedBlurLumi(x+1, y+1));
    gradY(x, y) =cast<float>( -clampedBlurLumi(x-1, y-1) + -clampedBlurLumi(x, y-1)*2
                              -clampedBlurLumi(x+1, y-1)
                             + clampedBlurLumi(x-1, y+1) + 2*clampedBlurLumi(x, y+1)
                             + clampedBlurLumi(x+1, y+1));

    //compute tensor and blur it
    tensor(x, y, c) = 0.f;
    tensor(x, y, 0) = gradX(x, y) * gradX(x, y);
    tensor(x, y, 1) = gradX(x, y) * gradY(x, y);
    tensor(x, y, 2) = gradY(x, y) * gradY(x, y);
    tensorClamped(x, y, c) = tensor(x_clamped, y_clamped, c);
    blurTensorX(x, y, c) = sum(tensorClamped(x+(rx-radius), y, c) * gaussKernel(rx));
    clampedBlurTensorX(x, y, c) = blurTensorX(x_clamped, y_clamped, c);
    blurTensor(x, y, c) = sum(clampedBlurTensorX(x, y+(rx-radius), c) * gaussKernel(rx));

     //get smallest eigenvector (eigenvector w/ smallest eigenvalue)
    Expr trace = blurTensor(x, y, 0) + blurTensor(x, y, 2);
    Expr det = blurTensor(x, y, 0)*blurTensor(x, y, 2) - blurTensor(x, y, 1)*blurTensor(x, y, 1);
    Expr eval1 = trace*.5f + pow((pow(trace, 2)*.25f - det), .5f);
    Expr eval2 = trace*.5f - pow((pow(trace, 2)*.25f - det), .5f);
    Expr vecx = select(eval1 < eval2, eval1-blurTensor(x, y, 2), eval2-blurTensor(x, y, 2));
    Expr vecy = blurTensor(x, y, 1);
    output(x, y) = Halide::atan2(vecy, vecx) + 3.1415926f;
    
    //schedule
    tensor.parallel(c);
    tensor.compute_root();
    blurTensorX.parallel(c);
    blurTensorX.compute_root();
    blurTensor.parallel(c);
    blurTensor.compute_root();
    output.compute_root();
    
    return output;
}

Func painterlyHalide(const Halide::Tools::Image<float> &im, const Halide::Tools::Image<float> &texture, int N, int size, float noise) {
    
    Var x("x"), y("y"), c("c");
    Var xo("xo"), xi("xi"), yo("yo"), yi("yi");
    Func output("output");
    Func brushLarge("brushLarge");
    Func brushSmall("brushSmall");
    Func textureClamped("textureClamped");
    Func scaledTexture("scaledTexture");
    Func scaledTextureFourth("scaledTextureFourth");
    Func lumi("lumi");
    Func clampedLumiX("clampedLumiX");
    Func blurLumiX("blurLumiX");
    Func clampedBlurLumiX("clampedBlurLumiX");
    Func blurLumi("blurLumi");
    Func highPass("highPass");
    Func energy("energy");
    Func clampedEnergyX("clampedEnergyX");
    Func blurEnergyX("blurEnergyX");
    Func clampedBlurEnergyX("clampedBlurEnergyX");
    Func blurEnergy("blurEnergy");
    Func sharpnessMap("sharpnessMap");
    
    Expr x_clamped = clamp(x, 0, im.width()-1);
    Expr y_clamped = clamp(y, 0, im.height()-1);
    float sigmaG = 1.0;
    float truncate = 3.0;
    int nAngles = 36;
    
    //compute Gaussians for blur
    int radius = sigmaG * truncate;
    int fwidth = 2 * radius + 1;
    int factorSigma = 4;
    
    Func GKernelUnNorm("GKernelUnnorm");
    Func GKernelSum   ("GKernelSum");
    Func GKernel      ("GKernel");
    RDom rx(0, fwidth);
    GKernelUnNorm(x) = exp(-((x-radius)*(x-radius))/(2.0f*sigmaG*sigmaG));
    GKernelSum   (x) = sum(GKernelUnNorm(rx));
    GKernel      (x) = GKernelUnNorm(x)/GKernelSum(0);
    
    GKernelUnNorm.compute_root();
    GKernelSum   .compute_root();
    GKernel      .compute_root();
    Halide::Tools::Image<float> kernel = GKernel.realize(fwidth);
    
    //kernel with 4sigma
    int radius2 = factorSigma * sigmaG * truncate;
    int fwidth2 = 2 * radius2 + 1;
    
    Func GKernelUnNorm2("GKernelUnnorm2");
    Func GKernelSum2   ("GKernelSum2");
    Func GKernel2      ("GKernel2");
    RDom rx2(0, fwidth2);
    GKernelUnNorm2(x) = exp(-((x-radius2)*(x-radius2))/(2.0f*sigmaG*sigmaG*factorSigma*factorSigma));
    GKernelSum2   (x) = sum(GKernelUnNorm2(rx2));
    GKernel2      (x) = GKernelUnNorm2(x)/GKernelSum2(0);
    
    GKernelUnNorm2.compute_root();
    GKernelSum2   .compute_root();
    GKernel2      .compute_root();
    Halide::Tools::Image<float> kernel2 = GKernel2.realize(fwidth2);
    
    //rescale texture
    int maxDim = max(texture.width(), texture.height());
    float factor = (float)(size) / maxDim;
    textureClamped(x, y) = texture(clamp(x, 0, texture.width()-1), clamp(y, 0, texture.height()-1));
    scaledTexture(x, y) = textureClamped(cast<int>(round(1/factor * x)),
                                         cast<int>(round(1/factor * y)));
    scaledTextureFourth(x, y) = textureClamped(cast<int>(round(4/factor * x)),
                                               cast<int>(round(4/factor * y)));
    
    int w = im.width(), h = im.height();
    int halfTexWidth = texture.width()*factor/2;
    int halfTexHeight = texture.height()*factor/2;
    RDom r(-halfTexWidth, texture.width()*factor, -halfTexHeight, texture.height()*factor, 0, N);
    
    //blurred luminance for sharpness mask
    lumi(x, y) = .3f*im(x, y, 0) + .6f*im(x, y, 1) + .1f*im(x, y, 2);
    clampedLumiX(x, y) = lumi(x_clamped, y_clamped);
    blurLumiX(x, y) = sum(clampedLumiX(x+(rx-radius), y) * GKernel(rx));
    clampedBlurLumiX(x, y) = blurLumiX(x_clamped, y_clamped);
    blurLumi(x, y) = sum(clampedBlurLumiX(x, y+(rx-radius)) * GKernel(rx));
    //blurred energy
    highPass(x, y) = lumi(x, y) - blurLumi(x, y);
    energy(x, y) = pow(highPass(x, y), 2);
    clampedEnergyX(x, y) = energy(x_clamped, y_clamped);
    blurEnergyX(x, y) = sum(clampedEnergyX(x+(rx2-radius2), y) * GKernel2(rx2));
    clampedBlurEnergyX(x, y) = blurEnergyX(x_clamped, y_clamped);
    blurEnergy(x, y) = sum(clampedBlurEnergyX(x, y+(rx2-radius2)) * GKernel2(rx2));
    
    //schedule importance map computation
    lumi.compute_root();
    blurLumiX.compute_root();
    blurLumi.compute_root();
    energy.compute_root();
    blurEnergyX.compute_root();
    blurEnergy.compute_root();
    RDom rIm(0, w, 0, h);
    Func blurEnergyMax;
    blurEnergyMax.compute_root();
    
    Expr energyDom = blurEnergy(rIm.x, rIm.y);
    blurEnergyMax(x) = maximum(energyDom);
    sharpnessMap(x, y) = blurEnergy(x, y)/blurEnergyMax(0);
    
    Halide::Tools::Image<float> importanceMap = sharpnessMap.realize(w, h);
    Expr importanceDom = importanceMap(rIm.x, rIm.y);
    int Nscaled = N/evaluate<float>(sum(importanceDom)/(float)(w*h));
    
    //rotate brushes
    Func angles = computeAnglesHalide(blurLumi, kernel2, w, h, radius2);
    Func rotatedBrushes = rotateBrushesHalide(scaledTexture, nAngles, halfTexWidth*2, halfTexHeight*2);
    Func rotatedBrushesFourth = rotateBrushesHalide(scaledTextureFourth, nAngles, halfTexWidth/2, halfTexHeight/2);

    //compute random locations for large and small detail
    Expr randX = cast<int> (random_float()*(w-halfTexWidth*2) + halfTexWidth);
    Expr randY = cast<int> (random_float()*(h-halfTexHeight*2) + halfTexHeight);
    Func noiseModLarge;
    noiseModLarge(x) = 1-noise/2 + noise*(random_float());
    
    Func randomLarge;
    randomLarge(x) = {randX, randY};
    
    brushLarge(x, y, c) = 0.f;
    Expr brushNumLarge = cast<int> (angles(randomLarge(r.z)[0], randomLarge(r.z)[1])*nAngles/(2*3.1415926f)) % nAngles;
    //r.z is range 0-N
    brushLarge(randomLarge(r.z)[0]+r.x, randomLarge(r.z)[1]+r.y, c) =
        lerp(brushLarge(randomLarge(r.z)[0]+r.x, randomLarge(r.z)[1]+r.y, c),
             im(randomLarge(r.z)[0], randomLarge(r.z)[1], c) * noiseModLarge(r.z),
             rotatedBrushes(r.x+halfTexWidth, r.y+halfTexHeight, brushNumLarge));
    
    //(1) compute Nscaled (x,y) pairs
    //(2) go through and randomly reject based on importance map, store these
    //(3) draw smaller brush at all of these locations
    RDom r2(-halfTexWidth/4, texture.width()*factor/4, -halfTexHeight/4, texture.height()*factor/4, 0, Nscaled);
    Func randomSmall, reject, noiseModSmall;
    randomSmall(x) = {randX, randY};
    reject(x) = random_float() < importanceMap(randomSmall(x)[0], randomSmall(x)[1]);
    noiseModSmall(x) = 1-noise/2 + noise*(random_float());
    
    brushSmall(x, y, c) = brushLarge(x, y, c);
    Expr brushNumSmall = cast<int> (angles(randomSmall(r2.z)[0], randomSmall(r2.z)[1])*nAngles/(2*3.1415926f)) % nAngles;
    Expr interpolatedVal =
        lerp(brushSmall(randomSmall(r2.z)[0]+r2.x, randomSmall(r2.z)[1]+r2.y, c),
             im(randomSmall(r2.z)[0], randomSmall(r2.z)[1], c) * noiseModSmall(r2.z),
             rotatedBrushesFourth(r2.x+halfTexWidth/4, r2.y+halfTexHeight/4, brushNumSmall));
    Expr noEffect = brushSmall(randomSmall(r2.z)[0]+r2.x, randomSmall(r2.z)[1]+r2.y, c);
    
    brushSmall(randomSmall(r2.z)[0]+r2.x, randomSmall(r2.z)[1]+r2.y, c) = select(reject(r2.z),
                                                                                 interpolatedVal, noEffect);

    scaledTexture.compute_root();
    scaledTextureFourth.compute_root();
    brushLarge.compute_root();
    
    return brushSmall;
}

