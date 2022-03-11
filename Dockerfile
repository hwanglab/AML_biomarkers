# syntax = docker/dockerfile:1.3
FROM rocker/r-ver:4.1.0

ENV RENV_VERSION 0.15.4-4

WORKDIR /src

RUN apt-get update
RUN apt-get install -y libcurl4-openssl-dev libhdf5-dev \
    imagemagick libfftw3-dev libgdal-dev proj-bin libgeos++-dev \
    libcairo2-dev libxt-dev libharfbuzz-dev libfribidi-dev \
    wget libgit2-dev default-jdk r-cran-rjava \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN --mount=type=secret,id=PAT \
    GITHUB_PAT=$(cat /run/secrets/PAT) \
    R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

ADD renv.lock renv.lock
ADD environment.yml environment.yml

RUN R CMD javareconf

ENV RENV_PATHS_CACHE /buildx_cache/renv
ENV RENV_CONFIG_CACHE_SYMLINKS FALSE

RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y python3.7
#RUN python3 -m venv python_venv 

#RUN bash setup_venv.sh
RUN mkdir renv
ADD renv renv

RUN cp renv/.Rprofile .Rprofile

RUN --mount=type=secret,id=PAT \
    --mount=type=cache,target=/buildx_cache/renv \
    GITHUB_PAT=$(cat /run/secrets/PAT) \
    R -e 'renv::activate();renv::restore();renv::isolate()'

RUN wget https://github.com/codercom/code-server/releases/download/1.408-vsc1.32.0/code-server1.408-vsc1.32.0-linux-x64.tar.gz \
     && tar -xzvf code-server1.408-vsc1.32.0-linux-x64.tar.gz \
     && chmod +x code-server1.408-vsc1.32.0-linux-x64/code-server

COPY cli cli
COPY lib lib
COPY docker_entry.sh docker_entry.sh

RUN mkdir clininfo
RUN mkdir outs
RUN mkdir plots

EXPOSE 80

ENTRYPOINT ["bash", "docker_entry.sh"]