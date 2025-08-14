# ABC-DS: obstacle Avoidance with Barrier-Certified polynomial Dynamical Systems

Accompanying code for the paper [Learning Barrier-Certified Polynomial Dynamical Systems for Obstacle Avoidance with Robots](https://ieeexplore.ieee.org/document/10610828) by
Martin Schonger<sup>1*</sup>,
Hugo T. M. Kussaba<sup>1*</sup>,
Lingyun Chen<sup>1</sup>,
Luis Figueredo<sup>2</sup>,
Abdalla Swikir<sup>1</sup>, 
Aude Billard<sup>3</sup>,
and Sami Haddadin<sup>1</sup>, published in the proceedings of the 2024 International Conference on Robotics and Automation (ICRA 2024). Pre-print available at [arXiv:2403.08178](https://arxiv.org/abs/2403.08178).

<sup>1</sup>Munich Institute of Robotics and Machine Intelligence (MIRMI), Technical University of Munich (TUM), Germany. Abdalla Swikir is also with the Department of Electrical and Electronic Engineering, Omar Al-Mukhtar University (OMU), Albaida, Libya.\
<sup>2</sup>School of Computer Science, University of Nottingham, UK. Luis Figueredo is also an Associated Fellow at the MIRMI, TUM.\
<sup>3</sup>Learning Algorithms and Systems Laboratory, EPFL, Switzerland.\
<sup>*</sup>These authors contributed equally to the paper.

Modified by Simone Silenzi to avoid PENBMI.

### Setup
Install MATLAB (tested with R2023a).

Install the required third party tools:
[YALMIP](https://yalmip.github.io/),
[PENLAB](https://web.mat.bham.ac.uk/kocvara/penlab/),
(Optional for plotting: [crameri colormaps](https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)).

> **Note**
> Make sure that the non-toolbox paths are before/on top of the toolbox paths.

Run:
```bash
git clone https://github.com/ssilenzi/abc-ds.git
cd abc-ds
git submodule init
git submodule update
```

### Usage
Open the `abc-ds` folder in MATLAB.

Configure the desired experiments in `main2.m` and run this script.

Check the `output` folder for results and logs.

(Optionally, recreate the plots from the paper with `generate_plots.m`, and the animations from the video with `generate_plots_video.m`.)


### Contact
martin.schonger@tum.de


This software was created as part of Martin Schonger's master's thesis in Computer Science at the Technical University of Munich's (TUM) School of Computation, Information and Technology (CIT).


Copyright © 2023 Martin Schonger
Copyright © 2025 Simone Silenzi
This software is licensed under the GPLv3.
