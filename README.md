# Granule-In-Cell (2D)
the official implementation of *The Granule-In-Cell Method for Simulating Sand–Water Mixtures*

![output](./results/bigball/output.gif)

Here we show a 2D example of sand and water simlation driven by Granule-In-Cell method. The engine is implemented by C++ built by xmake. As the total particle num in 2d scene is not so huge, the code is not fully parallelized. The viewer is written by taichi. To play with the code, you need to install one of the following compilers:

  - Windows: MSVC/MinGW
  - Linux: G++ (≥11.0) / Clang (≥13.0)
  - macOS: Apple Clang (≥13.0)

To build the engine, you need to install [xmake](https://xmake.io/) as instructed in https://xmake.io/#/getting_started.
Then, run the following commands successively in this directory for getting started:
```shell
xmake
xmake r demo -t falling -r 50 -e 201 -s 64
```
The exported images are subsequently generated in `build/[Platform]/[Arch]/release/output`.

To view the result, you need to install taichi based on python3.10:
```shell
pip install numpy imageio taichi
```
Then you can view your simultion by following commands:
```shell
python viewer.py -t falling -s 64
```
Finally you can get the image and the gif output in the `results/falling`

Also we provide an integrated operation method you can simulate and view by following commands:
```shell
python GIC.py
```


