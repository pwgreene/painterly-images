#include <iostream>

#include "a10.h"
#include "basicImageManipulation.h"
#include "timing.h"
#include <sstream>

using namespace std;


Image timeFunction(const Image &input, const Image &brush, Image (*fn)(const Image&, const Image&, int, int, float), int N, int size, int noise) {
    int w = input.width(), h = input.height();
    Image output(w, h, input.channels());
    unsigned long s = millisecond_timer();
    for (int i=0; i<N_TIMES; i++) {
        output = fn(input, brush, N, size, noise);
    }
    float total_time = float(millisecond_timer()-s);
    float mpixels = float(w*h)/1e6;
    cout << "Reference Implementation" << endl
    << "------------------------" << endl
    << "  - runtime " << total_time/N_TIMES << " ms " << endl
    << "  - throughput " << (mpixels*N_TIMES)/(total_time/1000) << " megapixels/sec" << endl;
    return output;
}

Halide::Tools::Image<float> timeFunction(const Halide::Tools::Image<float> &input,
                                         const Halide::Tools::Image<float> &brush,
                                         Func (*fn)(const Halide::Tools::Image<float>&, const Halide::Tools::Image<float>&, int, int, float),
                                         int N, int size, int noise) {
    int w = input.width(), h = input.height();
    Halide::Tools::Image<float> output;
    Func halideFunc;
    unsigned long s = millisecond_timer();
    for (int i=0; i<N_TIMES; i++) {
        halideFunc = fn(input, brush, N, size, noise);
        output = halideFunc.realize(w, h, 3);
    }
    float total_time = float(millisecond_timer()-s);
    float mpixels = float(w*h)/1e6;
    cout << "Halide Implementation" << endl
    << "------------------------" << endl
    << "  - runtime " << total_time/N_TIMES << " ms " << endl
    << "  - throughput " << (mpixels*N_TIMES)/(total_time/1000) << " megapixels/sec" << endl;
    return output;
}

//tests the singleScalePaint function on input image with given brush and writes the image to disk
void testSingleScale(const Image &input, const Image &brush, int N, int size, int noise) {
    Image output(input.width(), input.height(), input.channels());
    Image importance(input.width(), input.height(), 1);
    importance.set_color(1);
    singleScalePaint(input, output, importance, brush, N, size, noise);
    output.write("./Output/SingleScaleTest.png");
    cout << "image written to ./Output/SingleScaleTest.png\n";
}

//computes the angles for input and draws them using the brush texture
void testAngles(const Image &input, const Image &brush, int size) {
    Image angles = computeAngles(input);
    Image output = drawAngles(angles, brush, size);
    output.write("./Output/AnglesTest.png");
    cout << "image written to ./Output/AnglesTest.png\n";
}

//test orientedPaint function on input using brush texture
void testOrientedPaint(const Image &input, const Image &brush, int N, int size) {
    Image output = orientedPaint(input, brush, N, size);
    output.write("./Output/OrientedPaintTest.png");
    cout << "image written to ./Output/OrientedPaintTest.png\n";
}

//writes sharpness mask of input to disk
void testSharpnessMap(const Image &input) {
    Image output = sharpnessMap(input);
    output.write("./Output/SharpnessMaskTest.png");
    cout << "image written to ./Output/SharpnessMaskTest.png\n";
}

//runs the halide version on two-scale oriented paint algorithm and writes to disk
void testPainterlyHalide(const Halide::Tools::Image<float> &input, const Halide::Tools::Image<float> brush, int N, int size, float noise) {
    int w = input.width(), h = input.height();
    Func halideFunc = painterlyHalide(input, brush, N, size, noise);
    Halide::Tools::Image<float> output = halideFunc.realize(w, h, 3);
    Halide::Tools::Image<float> norm_out = normalize_image(output);
    Halide::Tools::save_image(norm_out, "./Output/HalidePainterlyTest.png");
    cout << "image written to ./Output/HalidePainterlyTest.png\n";
}

//times the computation time of two-scale oriented painterly rendering for both
//c++ and halide implementations
void testCompareTimes(const Image &input, const Halide::Tools::Image<float> &inputHalide,
                      const Image &brush, const Halide::Tools::Image<float> &brushHalide,
                      int N, int size, float noise) {
    Halide::Tools::Image<float> HalideOut = timeFunction(inputHalide, brushHalide, painterlyHalide, N, size, noise);
    Image refOut = timeFunction(input, brush, orientedPaint, N, size, noise);
}

//times various values of N on two of the same images
void testCompareValuesOfN(const Image &input, const Halide::Tools::Image<float> &inputHalide,
                          const Image &brush, const Halide::Tools::Image<float> &brushHalide,
                          int maxN, int size, float noise) {
    cout << "Testing with different values of N\n";
    for (int n = maxN/10; n <= maxN; n+= maxN/10) {
        cout << endl;
        cout << "testing with n = " << n << "\n";
        Halide::Tools::Image<float> HalideOut = timeFunction(inputHalide, brushHalide, painterlyHalide, n, size, noise);
        Image refOut = timeFunction(input, brush, orientedPaint, n, size, noise);
        //uncomment below to write to disk
//        std::ostringstream filenameRef
//        std::ostringstream filenameHalide
//        Halide::Tools::Image<float> norm_out = normalize_image(HalideOut);
//        filenameHalide << "./Output/timeTests/halide/test-halide" << n << ".png";
//        Halide::Tools::save_image(norm_out, filename.str());
//        filenameRef << "./Output/timeTests/ref/test-ref" << n << ".png";
//        refOut.write(filenameRef);
    }
}

int main()
{
    // Test your intermediate functions
    Image boston = Image("./Input/boston.png");
    Image london = Image("./Input/london.png");
    Image ville = Image("./Input/villeperdue.png");
    Image brush1 = Image("./Input/brush.png");
    Image brush2 = Image("./Input/longBrush.png");
    Image brush3 = Image("./Input/longBrush2.png");
    Image china("./Input/china.png");
    
    Halide::Tools::Image<float> bostonHalide = Halide::Tools::load_image("./Input/boston.png");
    Halide::Tools::Image<float> chinaHalide = Halide::Tools::load_image("./Input/china.png");
    Halide::Tools::Image<float> brush1Halide = Halide::Tools::load_image("./Input/brush.png");
    Halide::Tools::Image<float> brush2Halide = Halide::Tools::load_image("./Input/longBrush.png");
    Halide::Tools::Image<float> brush3Halide = Halide::Tools::load_image("./Input/longBrush2.png");
    Halide::Tools::Image<float> londonHalide = Halide::Tools::load_image("./Input/london.png");
//    Halide::Tools::Image<float> hkHalide = Halide::Tools::load_image("./Input/hk.png");
    
    
    int testN = 10000, testSize = 30;
    float testNoise = 0.3;
    
    testSingleScale(ville, brush1, testN, testSize, testNoise);
    testAngles(china, brush2, 10);
    testOrientedPaint(ville, brush1, testN, testSize);
    testSharpnessMap(boston);
    testPainterlyHalide(chinaHalide, brush3Halide, 20000, 50, testNoise);
    
    testCompareTimes(china, chinaHalide, brush1, brush1Halide, 10000, 30, 0.3);
    
    
    return EXIT_SUCCESS;
}
