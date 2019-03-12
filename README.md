# psraytracing
Raytracing algorithm for one-way travel-time calculation in passive seismic or microseismic applications. Written in Julia 1.1.0. 

This is a sister package providing same functions with the Python3 package [psraytrace](https://github.com/uqzzhao/psraytrace).

* shooting( ) - the core function of raytracing algorithm, actually a 2d algorithm
* psraytrace( ) - the core function to call shooting function and adapt it to 3D case. Unlike calculating two-way travel-time in reflection seismic raytracing program, only one-way travel-time for transmimitted P- or S-waves is calculated. You may modify this code to calculate the PS wave travel-time.

## Usage
You may refer to the Jupyter Notebook example code as psraytracing.ipynb included in this package for how to use the core function and may adapt to your case by modifying these examples. 

## License
MIT License.

## Authorship
psraytracing was written in 2019. The sole contributor is [Zhengguang Zhao](https://www.researchgate.net/profile/Zhengguang_Zhao2), who now works for DeepListen on part-time basis.
