# BMOBench
Welcome to **B**lack-Box **M**ulti-Objective **O**ptimization **B**enchmarking (**BMOB**) Platform.

The aim of this platform is to /consolidate/ black-box multi-objectives problems from the literature into a single framework; which makes it easier for researchers in the Multi-Objective Optimization community to compare, assess, and analyze previous and new algorithms comprehensively. In essence, adding a brick to tools for *reproducible research*.

With *BMOBench*, you can test your newly developed algorithms on 100 established problems from the multi-objective optimization community and automatically get the experiments results in latex-based paper template. The results--as data profiles--are reported in terms of four quality indicators: hypervolume, additive epsilon-indicator, inverted generational distance, and generational distance.

# Requirments:

* For experiments: `MATLAB`, `C` [compiler](http://mingw-w64.org/doku.php) 
* For post-processing: `C` [compiler](http://mingw-w64.org/doku.php) and `Python 32-bit` with the Numpy, matplotlib, and palettable packages, [Anaconda](https://www.continuum.io/downloads) is a good start.


# Setup

At this point of time, the platform is supporting `MATLAB` and `C`. Future releases may support `Python` as well.

To start with the *BMOBench* platform, download its code from `github`:
* Download as [ZIP](https://github.com/ash-aldujaili/BMOBench/zipball/master)/ [TAR](https://github.com/ash-aldujaili/BMOBench/tarball/master)
* Unzip(tar) the downloaded file to find the following folders: 
  * `problems`: problems-related files and descriptions 
  * `postprocess`: for data post-processing 
  * `matlab` : for running experiments in `MATLAB` 
  * `c` : for running experiments in `C` 
  * `latex-template`: paper template incorporating results generated from the post-processing step.

### Experiments Setup

#### MATLAB testbed
* Before you can test your algorithm, you need to run this code for once, to compile some necessary `mex` files to speed up the computation:
1. From the MATLAB command window, `cd` to `matlab/benchmark`
2. Execute the `setup.m` script. 
~~~
>>setup
~~~

#### C testbed
* Before you can test your algorithm in C:
1. Put your algorithm in `c/algs` similar to the random search baseline algorithm, `MO-RANDOM`. 
2. Edit the self-explanatory `c/benchmark/main.c` file to include your algorithm. You may want to change the number of runs  (`NUM_RUNS`) as well as the budget multiplier factor (`BUDGET_MULTIPLIER`) for your experiments. This is can be done by editing the `c/benchmark/globaldeclare.h` file.



### Post-Processing Setup

* Before you can post-process your results, you need to run this code for once, to compile some necessary `C` libraries to speed up the post-processing computation:
1. From the system shell, `cd` to `postprocess/scripts`
2. Execute the `setup.py` python script. 
~~~
>>python setup.py
~~~

# Experiments

## Running Experiments (MATLAB)
1. Launch MATLAB and set your current path to the `matlab` directory 
2. Put your favourite/developed algorithm under the directory `matlab/algs` directory
3. Incorporate the function call to your algorithms in the self-explanatory files `matlab/benchmark/run_performance_experiments.m` and `matlab/benchmark/run_timing_performance.m` in a way similar to the exemplar algorithm.
4. Edit the main script `matlab/runBenchmark.m` by incorporating the algorithm's name; If you wish to change the number of function evaluations and the number of runs, you can do this by editing the same.
5. Now, you are ready to run the experiments by simply running `matlab/runBenchmark.m`. The collected results will be generated into a new directory `EXP_RESULTS`.
6. Since computing the hypervolume (hv) indicator values is computationally expensive, we have made use of the readily available mex routines for that purpose. Edit `matlab/run_hv_computation.m` to include the algorithms of interest (whose results are now in the `EXP_RESULTS` directory). The computed hv value will also stored in the `EXP_RESULTS` directory.

## Running Experiments (C)

1. `cd` to the `c` directory in a system shell, and hit `make` to compile.
2. Execute the `runme` executable. The collected results will be generated into a new directory `EXP_RESULTS`.
3. You may want to do `step 6` of MATLAB Experiment procedure (above) to compute the hv profile.


## Post-Processing Data

If you have already ran the `postprocess/scripts/setup.py`, you can generate the results as follow:

1. From the system shell, `cd` to `postprocess/scripts`.
2. Edit the `data` dictionary within the main scripting file `run_postprocessing` to incorporate the benchmarked algorithms similar to the description of the exemplar algorithm, you may choose among the available quality indicators to compare. 
  * All the quality indicators are computed in this stage except for the hypervolume (hv) indicator. If you want to include it in the processing, you have to execute `step 6` of MATLAB Experiment procedure (above). 
  * The platform also allows comparing deterministic vs stochastic algorithms by reporting the best data profiles from the multi-run (stochastic) algorithms rather than their means against the single-run data profile of the deterministic algorithm which is usually given the same evaluation budget/run times the number of runs given to the stochastic algorithms. To report the mean data profiles set `isMean` to *True*. To report the best data profiles, set it to *False*. For more details, please refer to the technical report.
3. Execute the `run_postprocessing.py` python script.  
~~~
>>python run_postprocessing.py
~~~


## Compiling the LaTeX paper.

With a LaTeX editor, the paper can be compiled directly. Note than the paper only compiles a portion of the results: *timing* and the *aggregated* performance of the algorithms over the problem categories and used quality indicators. The rest of the generated results can be found in `postprocess/postproc`.


# Acknowledgement / Inspiration
This platform is inspired by two solid papers from the literature:

1. Custódio, Ana Luísa, et al. "**[Direct multisearch for multiobjective optimization.](http://www.mat.uc.pt/~lnv/papers/dms.pdf)**" *SIAM Journal on Optimization* 21.3 (2011): 1109-1140.

2. Brockhoff, Dimo, Thanh-Do Tran, and Nikolaus Hansen. "**[Benchmarking numerical multiobjective optimizers revisited.](https://hal.inria.fr/hal-01146741/document)**" *Genetic and Evolutionary Computation Conference (GECCO 2015)*. 2015.

**Thanks** extended to *[Bhavarth Pandya](https://github.com/bhavarthpandya), [Chaitanya Prasad](https://github.com/chaitanya94), [Khyati Mahajan](https://github.com/khyatimahajan), and [Shaleen Gupta](https://github.com/shaleenx)*. They contributed greatly to the core of the `C` platform and verified the correctness of the `C`-coded problems to their `AMPL` counterpart.

# Citation

If you write a scientific paper describing research that made use of this code, please cite the following [report](http://arxiv.org/abs/1605.07009):

~~~
@article{bmobench-16,
  author = {Abdullah Al-Dujaili and S. Suresh},
  title = {BMOBench: Black-Box Multi-Objective Optimization Benchmarking Platform},
  journal = {ArXiv e-prints},
  year = {2016},
  volume = {arXiv:1605.07009},
  url = {http://arxiv.org/abs/1605.07009}
}
~~~
