FROM rocker/r-ver:4.5.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    unzip \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libtbb-dev \
    libblas-dev \
    liblapack-dev \
    libxt-dev \
    r-base-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/FinemapColoc

WORKDIR /opt/FinemapColoc

RUN Rscript deps/install.R

RUN curl -fsSL https://install.julialang.org | sh -s -- --yes

ENV PATH="/root/.juliaup/bin:${PATH}"

RUN juliaup add 1.11 && juliaup default 1.11

RUN R_HOME=/usr/local/lib/R julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
