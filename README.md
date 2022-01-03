# Benchmark UNU.RAN C and Python API

Benchmarks for the paper "Automatic random variate generation in Python"

### Build Instructions

Build UNU.RAN and Rngstreams (optional, build Rngstreams only if you want to use it) and change the paths in Makefile according to your installation.
By default, the default expected installation path is `./unuran-1.8.1/build` and `./rngstreams-1.0.1/build` respectively.

After building UNU.RAN and Rngstreams (optionally), run `make all` or `make` to build the C benchmark code and execute it.

To run Python benchmarks, execute the `vp_benchmark.ipynb` notebook in Jupyter. Following dependencies must be installed to run python benchmarks:

- `scipy==1.8.1`
- `numpy>1.19.0`
