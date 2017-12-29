# painterly-images

Creates non-photorealistic renderings of images using the Halide image processing library for fast, scalable compuatations. 

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

