
# C++11 Ray Tracing / Marching Framework

This is a basic C++11 OpenGL viewer application running a ray tracing / marching kernel for each tile in a pool of worker threads. The purpose of this application is to be a test framework for some rendering ideas I plan to explore.

# Features

A basic rendering framework and a collection of different libraries.

* Point, ray and triangle with AABB intersection
* Two different ray / triangle intersection tests
* Triangle / line segment distance functions
* Ray / plane intersection
* Template classes for vectors and matrices with common computer graphics operations implemented
* OpenGL viewer with text rendering and frame limiter
* Ten runtime switchable test scenes
* Screenshot feature, writing of BMP images
* Ray generation for orthographic and perspective cameras
* Build-in Cornell Box scene
* Tiled, parallel, CPU writeable + OpenGL drawable frame buffer system
* 3D DDA based grid intersection accelerator
* Basic mesh processing pipeline (load from disk, transform, compute normals etc.)
* Sampling module supporting various Low Discrepancy Sequences / QMC methods
* Time helpers for logging and performance measurement
* Basic tracing / logging system 
* Super-sampling
* Gamma correction

# Building

Compiling has no prerequisites aside from a C++11 ready Clang and GLUT. I tested with Apple's [command line development tools][devdownloads]. Building should succeed with a simple

[devdownloads]:https://developer.apple.com/downloads/

    $ make

# Images

Not much to see, just the basic are in place.

![viewer](https://raw.github.com/blitzcode/cpp-11-ray-trace-march-framework/master/img/viewer.png)
![scenes](https://raw.github.com/blitzcode/cpp-11-ray-trace-march-framework/master/img/scenes.png)

# Legal

This program is published under the [MIT License](http://en.wikipedia.org/wiki/MIT_License).

# Author

Developed by Tim C. Schroeder, visit my [website](http://www.blitzcode.net) to learn more.

