# RECAP: Detecting Evolutionary Patterns of Cancers using Consensus Trees

RECAP is an iterative, heuristic algorithm for solving the Mulitple Choice Consensus Tree (MCCT) problem. The input to RECAP is a family of sets of patient tumor phylogenies and an integer k > 0. RECAP then returns (i) a single expanded tumor phylogeny for each patient, (ii) an assignment of patients into k clusters, and (iii) a consensus tree for each cluster summarizing its repeated evolutionary patterns.  

![RECAP solves the Multiple Choice Consensus Tree problem.](./overview.png)

RECAP can be run over a range of values for the number of clusters, k, and then a model selection criterion, such as the elbow method or Bayesian Information Criterion, can be used to select a particular k. 

## Contents

  1. [Getting started](#start)
     * [Dependencies](#dep)
     * [Compiling code](#build)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
     * [Running RECAP](#recap)
     * [Running Analyze](#analyze)

<a name="start"></a>

## Getting started

RECAP is written in C++11. The repository is organized as follows:

| Folder    | Description                                                  |
| --------- | ------------------------------------------------------------ |
| `src`     | source code for RECAP                                        |
| `data`    | input data for RECAP                                         |
| `results` | output from RECAP and example Python code to perform model selection |

<a name="dep"></a>

### Dependencies   

RECAP can be compiled with [CMake](https://cmake.org/) (>= 2.8) and has the following dependencies:

* [Boost](https://www.boost.org/) (>= 1.70.0)
* [Lemon](https://lemon.cs.elte.hu/trac/lemon) (>= 1.3.1)

<a name="build"></a>

## Compiling code

Here we walk through how to build the executables using CMake. First, navigate to a directory in which you would like to download RECAP. Then perform the following steps: 

```bash
# Download repository
git clone https://github.com/elkebir-group/RECAP.git

# Enter downloaded RECAP folder
cd RECAP/

# Make new build directory and enter it
mkdir build
cd ./build/

# Use cmake to compile executables. Specify lemon path if not detected automatically.
cmake .. -DLIBLEMON_ROOT=/usr/local/
make
```

<a name="usage"></a>

## Usage instructions

Here we describe how to run the RECAP executables. 

<a name="io"></a>

### I/O format

RECAP takes as input a file containing a cohort of patient trees formatted according to a specific convention. An example of such a file can be found in the data folder [here](./data/TRACERx_lung/tracerx_lung.txt). 

Note that **all patient trees must have the same root node**. The label used in this example is "GL", but this must be indicated on input to RECAP. 

Also, **mutation cluster labels should separate mutations with a consistent symbol**. The default is a semicolon ";" but this again can be indicated on input to RECAP. 

Finally, all **text following the pound symbol on each line will be ignored** by RECAP so it can be annotated in way that is useful to the user.      

The first line of the file should give the number of patients in the file.

```tex
99 # patients
```

On the second line, we indicate the number of trees for the first patient. 

```tex
3 #trees for CRUK0001
```

On the third line, we indicate the number of edges in the first tree for the first patient. 

```tex
4 #edges
```

The following lines each contain one directed edge from this first tree. The two node labels should be separated by a space. The first label is the source and the second label is the target for the directed edge. In this example, some of the node labels correspond to mutations clusters comprised of several mutations separated by a semicolon. 

```tex
GL TP53;MGA;WRN;EGFR
PASK ARHGAP35
TP53;MGA;WRN;EGFR PASK
ARHGAP35 NF1
```

The remaining trees for the first patient are indicated in a similar manner on the subsequent lines. 

```tex
4 #edges
GL TP53;MGA;WRN;EGFR
PASK ARHGAP35
TP53;MGA;WRN;EGFR NF1
NF1 PASK
4 #edges
GL TP53;MGA;WRN;EGFR
PASK ARHGAP35
PASK NF1
TP53;MGA;WRN;EGFR PASK
```

This format then repeats for the remaining 98 patients.  

```tex
1 #trees for CRUK0003
2 #edges
GL PIK3CA;EGFR;CDKN2A
PIK3CA;EGFR;CDKN2A CTNNB1
1 #trees for CRUK0004
2 #edges
GL TP53;EGFR
TP53;EGFR SMAD4;NOTCH1
1 #trees for CRUK0005
2 #edges
GL CMTR2;TP53;BRAF;PASK;TERT
CMTR2;TP53;BRAF;PASK;TERT NRAS
3 #trees for CRUK0006
3 #edges
GL PLXNB2;TP53;KEAP1;TERT
FANCC MAP3K1
PLXNB2;TP53;KEAP1;TERT FANCC
3 #edges
GL PLXNB2;TP53;KEAP1;TERT
MAP3K1 FANCC
PLXNB2;TP53;KEAP1;TERT MAP3K1
3 #edges
GL PLXNB2;TP53;KEAP1;TERT
PLXNB2;TP53;KEAP1;TERT FANCC
PLXNB2;TP53;KEAP1;TERT MAP3K1
.
.
.
1 #trees for CRUK0099
1 #edges
GL STK11;KEAP1;TP53
```

<a name="recap"></a>

### Running RECAP

This executable requires as input a filepath to the file of patient trees as described above. There are also several flags that must be set. The flag `-k`  should be followed with the number of desired clusters. The flag `-p` should be followed with a unique path for the output files. The flag `-R` should be followed with the label of the root node shared by all patient trees.   

```bash
# An example of the minimum input information needed to run RECAP
./build/recap -k 3 -p "./results/k3" -R "GL" "./data/TRACERx_lung/tracerx_lung.txt"
```

There are additional input options that can be specified.

`-r` allows the user to specify the number of restarts (default: 50). 

`-t` specifies a time limit in seconds as the stopping criterion (default: -1, no time limit). 

`-s` specifies a seed for the random number generator (default: 0). 

`-S` specifies the the separator used between mutations in a mutation cluster label (default: ;).  

`-c` specifies the minimum number of patients that must have a mutation for it to be included in the analysis. All mutations falling below this threshold are removed from the patient input trees (default: 1).  

These options can also be viewed by typing `./build/recap --help`.  

<a name="analyze"></a>

### Running Analyze

This executable takes as input the solution.txt file produced as output of RECAP. It also takes as input the original input file to RECAP described above in [I/O formats](#io). `-H` is an optional flag that indicates the column labels should be printed in the output.  

```bash
# Run analyze command on real data. 
./build/analyze "./results/TRACERx_lung/k3.solution.txt" \
-i "./data/TRACERx_lung/tracerx_lung.txt"
```

If the input instance was over simulated data, it is also possible to compare the output of RECAP to ground truth. This requires also indicating a ground truth solution file using `-s`. 

```bash
# Run analyze command on simulated data. 
./build/analyze "./results/simulations/RECAP/M5_m5/s0_k1_n100_M5_m5.solution.txt" \
-s "./data/simulations/M5_m5/simulations_solution/s0_k1_n100_M5_m5_simulations.txt" \
-i "./data/simulations/M5_m5/simulations_input/s0_k1_n100_M5_m5_simulations.txt"
```

In practice, it is often helpful to loop over all k values run for the input instance. This can be done with a bash script similar to the following:

```bash
#!/bin/bash
first=true 
for k in {1..15}
do
	f=./results/k${k}.solution.txt
	if [ "$first" = true ]
  	then
    	first=false
    	./build/analyze -H $f -i "./data/TRACERx_lung/tracerx_lung.txt"
    else
    	./build/analyze $f -i "./data/TRACERx_lung/tracerx_lung.txt"
    fi
done
```

The output of analyze can be saved to a file by redirecting it  `AnalyzeCommand > results.tsv `. 