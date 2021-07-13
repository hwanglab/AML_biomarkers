FROM rocker/r-ver:4.1.0

# set up system enviorment
RUN apt-get update
RUN apt-get install -y python3 python3-pip 
RUN apt-get install -y git
RUN apt-get install -y libcurl4-openssl-dev libssl-dev wget
RUN apt-get install -y pkg-config build-essential libgfortran-7-dev
RUN apt-get install -y libhdf5-dev libxml2-dev
RUN apt-get install -y x11-common libpng-dev libcairo2-dev libxt-dev
RUN apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN apt-get install -y libudunits2-dev libproj-dev
RUN apt-get install -y libharfbuzz-dev libfribidi-dev
RUN apt-get install -y xvfb
RUN apt-get install -y libgdal-dev tk
RUN apt-get install -y libmagick++-dev
RUN apt-get install -y fftw3-dev
RUN apt-get install -y default-jdk

# set up R enviroment
COPY renv.lock renv.lock
COPY /renv/activate.R /renv/activate.R
COPY requirements.txt requirements.txt

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir radian
#RUN Rscript -e "source('renv/activate.R');renv::restore(prompt = FALSE)"

# copy source code
COPY /lib /lib
COPY /cli /cli

# set up env vars
ENV AML_DATA=/sc_data

# set up folder structure
RUN mkdir -p /outs
RUN mkdir -p /plots
RUN mkdir -p /sc_data
RUN mkdir -p ~/.local/share/renv

RUN mkdir /data
RUN mkdir /clinical_info
