#!/bin/bash
#SBATCH --account=PCCF0022
#SBATCH --job-name=install_packages
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --output=slurm/latest_out.txt
#SBATCH --error=slurm/latest_err.txt

cd $SLURM_SUBMIT_DIR


LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/gdal/3.2.2/lib:/apps/proj/7.2.1/lib:/apps/geos/3.8.2/lib
PKG_CONFIG_PATH=/apps/fftw3/intel/19.0/mvapich2/2.3/3.3.8/lib/pkgconfig
GDAL_DATA=/apps/gdal/3.2.2/share/gdal

Rscript -e 'install.packages("devtools")'

Rscript -e 'install.packages("fftw", configure.vars = "PKG_CONFIG_PATH=/apps/fftw3/intel/19.0/mvapich2/2.3/3.3.8/lib/pkgconfig")'
#Rscript -e 'devtools::install_version("sf", version = "0.9-8", configure.args=c("--with-gdal-config=/apps/gdal/3.2.2/bin/gdal-config","--with-proj-include=/apps/proj/7.2.1/include","--with-proj-lib=/apps/proj/7.2.1/lib","--with-geos-config=/apps/geos/3.8.2/bin/geos-config"))'

Rscript -e 'dyn.load("/apps/gdal/3.2.2/lib/libgdal.so");dyn.load("/apps/geos/3.8.2/lib/libgeos_c.so", local=FALSE);devtools::install_version("spdep", version = "1.1-8");renv::restore()'
