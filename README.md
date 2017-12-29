# Halide Painterly Images

![reference image](london-reference.png)

Creates non-photorealistic renderings of images using the [Halide](halide-lang.org) image processing language for fast, scalable compuatations. 

To run from command line, cd to this directory and simply type

```
./painterly -i [INPUT PNG] -b [BRUSH PNG] -o [OUTPUT PNG]
```

replacing [INPUT PNG] with the image you would like to painterly render, [BRUSH PNG] with the brush image (some sample brushes are given in the Input/ directory), and [OUTPUT PNG] with the desired location of the output image.

Type 

```
./painterly --help
```

for the full usage.

If you would like to make modifications, you must have the [Halide](halide-lang.org) library on your machine, and modify the MakeFile to point HALIDE_LIB to the correct location.

# Why Halide?

The Halide implementation of this algorithm offers a much better performing computation than the vanilla C++, and is highly scalable, both in the size of the input image and the number of brush strokes. I spent a lot of time fine-tuning the scheduling, but there are still some improvements that could be done.

The following are results from running a comparison of the Halide implementation vs. standard C++ implementation on a 2014 2.6 GHz quad-core Macbook pro: 

| N       | 5000| 10000  | 15000  |20000  |30000  |40000  |50000  |
| -----|-----|-----|-----|-----|-----|-----|-----|
| Halide time (s)  | 2.128 | 2.202 |2.942 |3.598 |6.219 |6.055 |5.468 |
| C++ time (s)     | 27.67      |50.433| 60+ |60+ |60+ |60+ |60+ |

(Note: timing was cut off at 60 seconds)
