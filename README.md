# BMOBench
Welcome to **B**lack-Box **M**ulti-Objective **O**ptimization **B**enchmarking (**BMOB**) Platform.

The aim of this platform is to /consolidate/ black-box multi-objectives problems from the literature into a single framework; which makes it easier for researchers in the Multi-Objective Optimization community to compare, assess, and analyze previous and new algorithms comprehensively. In essence, adding a brick to tools for *reproducible research*.

With *BMOBench*, you can test your newly developed algorithms on 100 established problems from the multi-objective optimization community and automatically get the experiments results in latex-based paper template. The results--as data profiles--are reported in terms of four quality indicators: hypervolume, additive epsilon-indicator, inverted generational distance, and generational distance.

# Requirments:

* For experiments: `MATLAB`
* For post-processing: `C` [compiler](http://mingw-w64.org/doku.php) and `Python 32-bit` with the Numpy, matplotlib, and palettable packages, [Anaconda](https://www.continuum.io/downloads) is a good start.


# Setup

At this point of time, the platform is only supporting `MATLAB`. We are currently working on getting the `C` version up and running as soon as possible.

To start with the *BMOBench* platform, download its code from `github`:
* Download as [ZIP](https://github.com/ash-aldujaili/BMOBench/zipball/master)/ [TAR](https://github.com/ash-aldujaili/BMOBench/tarball/master)
* Unzip(tar) the downloaded file to find the following folders: 
  * `problems`: problems-related files and descriptions 
  * `postprocess`: for data post-processing 
  * `matlab` : for running experiments in `MATLAB` 
  * `latex-template`: paper template incorporating results generated from the post-processing step.

### Experiments Setup

* Before you can test your algorithm, you need to run this code for once, to compile some necessary `mex` files to speed up the computation:
1. From the MATLAB command window, `cd` to `matlab/benchmark`
2. Execute the `setup.m` script. 
~~~
>>setup
~~~

### Post-Processing Setup

* Before you can post-process your results, you need to run this code for once, to compile some necessary `C` libraries to speed up the post-processing computation:
1. From the system shell, `cd` to `postprocess/scripts`
2. Execute the `setup.py` python script. 
~~~
>>python setup.py
~~~

## Running Experiments (MATLAB)
1. Launch MATLAB and set your current path to the `matlab` directory 
2. Put your favourite/developed algorithm under the directory `matlab/algs` directory
3. Incorporate the function call to your algorithms in the self-explanatory files `matlab/benchmark/run_performance_experiments.m` and `matlab/benchmark/run_timing_performance.m` in a way similar to the exemplar algorithm.
4. Edit the main script `matlab/runBenchmark.m` by incorporating the algorithm's name; If you wish to change the number of function evaluations and the number of runs, you can do this by editing the same.
5. Now, you are ready to run the experiments by simply running `matlab/runBenchmark.m`. The collected results will be generated into a new directory `EXP_RESULTS`.





## Post-Processing Data

If you have already ran the `postprocess/scripts/setup.py`, you can generate the results as follow:
1. From the system shell, `cd` to `postprocess/scripts`.

2. Edit the `data` dictionary within the main scripting file `run_postprocessing` to incorporate the benchmarked algorithms similar to the description of the exemplar algorithm, you may choose among the available quality indicators to compare. 
  * The platform also allows comparing deterministic vs stochastic algorithms by reporting the best data profiles from the multi-run (stochastic) algorithms rather than their means against the single-run data profile of the deterministic algorithm which is usually given the same evaluation budget/run times the number of runs given to the stochastic algorithms. To report the mean data profiles set `isMean` to *True*. To report the best data profiles, set it to *False*. For more details, please refer to the technical report.

3. Execute the `run_postprocessing.py` python script.  
~~~
>>python run_postprocessing.py
~~~


## Compiling the LaTeX paper.

With a LaTeX editor, the paper can be compiled directly. Note than the paper only compiles a portion of the results: *timing* and the *aggregated* performance of the algorithms over the problem categories and used quality indicators. The rest of the generated results can be found in `postprocess/postproc`.


# Inspiration
This platform is inspired by two solid papers from the literature:

1. Custódio, Ana Luísa, et al. "**[Direct multisearch for multiobjective optimization.](http://www.mat.uc.pt/~lnv/papers/dms.pdf)**" *SIAM Journal on Optimization* 21.3 (2011): 1109-1140.

2. Brockhoff, Dimo, Thanh-Do Tran, and Nikolaus Hansen. "**[Benchmarking numerical multiobjective optimizers revisited.](https://hal.inria.fr/hal-01146741/document)**" *Genetic and Evolutionary Computation Conference (GECCO 2015)*. 2015.
