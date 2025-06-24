# MS-Thesis

This repository contains the following information for my Master's thesis "Effects of electron emission from biased electrodes on sheath dynamics under fusion-relevant conditions":

- Scripts for calculating and plotting the ion- and electron-induced electrom emission yield and spectra
- Experimental data for the ion- and electron-induced yield and spectra
- Gkeyll input files for simulating plasma sheaths with ion- and electron-induced emission at various bias potentials
- Outputs from Gkeyll for the simulations presented in the thesis document
- Postprocessing scripts to generate plots from .gkyl output files

This project uses the Gkeyll software (https://gkeyll.readthedocs.io/en/latest/index.html) to simulate the plasma sheath near electron emitting surfaces. In particular these surfaces are electrodes used to drive current through the plasma by applying a large voltage difference across them. This framework is applied to study the sheared-flow Z-pinch fusion concept. To reproduce the results presented in the thesis document, one must first install Gkeyll on a local or remote machine. The version of Gkeyll used here is called gkylzero, which can be pulled from this repository: https://github.com/ammarhakim/gkylzero/tree/main. To install gkylzero as well as the postprocessing tool postgkyl, do the following:

GKYLZERO
- From a linux terminal (WSL on windows for example) clone the repository with, git clone https://github.com/ammarhakim/gkylzero.git
- Navigate to gkylzero directory, cd gkylzero
- Open the file mkdeps.[SYSTEM].sh with, nano ./machines/mkdeps.[SYSTEM].sh. The system is whatever operating system you're running on, so linux for WSL and macos for mac.
- On the same line as the other flags (--build-openblas=yes --build-superlu=yes) add --build-openmpi=yes. This tells the computer to install Open MPI along with the other dependencies.
- Exit and run the file with, ./machines/mkdeps.[SYSTEM].sh, to install the dependencies.
- Open the file configure.[SYSTEM].cpu.sh with, nano ./machines/configure.[SYSTEM].cpu.sh
- Add the following to your file on the third line: --use-mpi=yes
- Make sure you're in the gkylzero directory and then run ./machines/configure.[SYSTEM].cpu.sh. In the output you should see USE_MPI=1 as well as filepaths to the include and library directories of Open MPI.
- Now install Gekyll with make install -j #. Replace the # with the number of cores you'd like to install with.

POSTGKYL
- Download miniconda with, wget https://repo.anaconda.com/miniconda/Miniconda3-py311_24.4.0-0-Linux-x86_64.sh -O /path/to/miniconda3/miniconda.sh
- Run, bash /path/to/miniconda3/miniconda.sh -b -u -p /path/to/miniconda3 and rm -rf /path/to/miniconda3/miniconda.sh
- Initialize with, /path/to/miniconda3/bin/conda init bash and /path/to/miniconda3/bin/conda init zsh
- Now clone the repository with, git clone https://github.com/ammarhakim/postgkyl.git
- Navigate to postgkyl directory, cd postgkyl
- Run the following: conda env create -f environment.yml, conda activate pgkyl, and pip install -e .
- You can activate postgkyl after reloading your terminal or doing source .bashrc with, conda activate pgkyl (this will need to be done every time before using postgkyl after loading up your terminal unless this command is added directly to your .bashrc file)
- Postgkyl also installs a built-in python library that will be used in the postprocessing scripts here
