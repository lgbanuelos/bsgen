# BsGen: A generalization measure for discovered process models based on bootstrapping

Process mining studies ways to derive value from process executions recorded in event logs of IT-systems, with
process discovery the task of inferring a process model for an event log emitted by some unknown system. One quality criterion
for discovered process models is generalization. Generalization seeks to quantify how well the discovered model describes future
executions of the system, and is perhaps the least understood quality criterion in process mining. The lack of understanding
is primarily a consequence of generalization seeking to measure properties over the entire future behavior of the system, when the
only available sample of behavior is that provided by the event log itself. In this paper, we draw inspiration from computational
statistics, and employ a bootstrap approach to estimate properties of a population based on a sample. Specifically, we define an
estimator of the model’s generalization based on the event log it was discovered from, and then use bootstrapping to measure the
generalization of the model with respect to the system, and its statistical significance. Experiments demonstrate the feasibility
of the approach in industrial settings.

This Git repository contains the code, implementing the method and experiment described in the paper:

Artem Polyvyanyy, Alistair Moffat and Luciano García-Bañuelos, Bootstrapping Generalization of Process Models Discovered From Event Data, [CoRR abs/2107.03876](https://arxiv.org/abs/2107.03876) (2021)

## Software requirements

- JDK 1.8  (Preferrably Oracle's JDK, as we found variations in the results when using the OpenJDK)
- Python 3
- PM4Py

## Usage

Clone this repository, e.g.

```
git clone https://github.com/lgbanuelos/bsgen
```

Given the large number of parameters to be set, it is preferrable that you set them directly in the script. To that end, take a look the lines 224 and onwards and check the comments. As you will see, it is assume the input Petri nets should be located in a folder which must be specified in line 225 (variable `input_dir`). Provide also a target directory setting the variable `output_dir` accordingly. Assumming you have all the required software properly installed, you would just need to execute from a terminal window the following command:

```
python3 bsgen.py
```

Note that the computation is long lasting. The experiment reported in the paper run for about 40 hours in a computer with 32 cores (the python script will allocate as many threads as cores available).

