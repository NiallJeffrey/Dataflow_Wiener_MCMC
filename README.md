# Dataflow_Wiener_MCMC

We use Maxeler Dataflow Engines (DFE) to construct an efficient Wiener filter of noisy and incomplete image data, and to quickly draw probabilistic samples of the compatible true underlying images from the Wiener posterior. Dataflow computing is a powerful approach using reconfigurable hardware, which can be deeply pipelined and is intrinsically parallel.

The Wiener-filtered image is the minimum-variance linear estimate of the true image (if the signal and noise covariances are known) and the most probable true image (if the signal and noise are Gaussian distributed). Many images are compatible with the data with different probabilities, given by the analytic posterior probability distribution referred to as the Wiener posterior.} The DFE code draws large numbers of samples of true images from this posterior. Naive computation of the Wiener-filtered image is impractical for large datasets, as it scales as $n^3$, where $n$ is the number of pixels. We use a messenger field algorithm, which is well suited to a DFE implementation, to draw samples from the Wiener posterior, that is, with the correct probability we draw samples of noiseless images that are compatible with the observed noisy image. 

## How to compile on MPC-X (Hartree specific)

```bash
#!/bin/sh
# name
#$ -N lr_compile
# shell
#$ -S /bin/sh
# execute script from the current directory
#$ -cwd
# set up the necessary environment
source /panfs/panfs.maxeler/maxeler/maxcompiler-2014.1.1/settings.sh
export PATH=/panfs/panfs.maxeler/maxeler/altera/13.1/quartus/bin:${PATH}
export MAXELEROSDIR=/opt/maxeler/maxeleros
# build the code
pushd {$repository}/RunRules/DFE
# set display to :99 because that is where Xvfb is running
export DISPLAY=":99"
make RUNRULE=DFE build
popd
```
