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

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.9-linux-x86_64.tar.gz \
    && tar -xzf julia-1.11.9-linux-x86_64.tar.gz \
    && rm julia-1.11.9-linux-x86_64.tar.gz \
    && mv julia-1.11.9 /opt

ENV PATH="/opt/julia-1.11.9/bin/:${PATH}"

ENV JULIA_DEPOT_PATH=/opt/julia_depot

RUN R_HOME=/usr/local/lib/R julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
