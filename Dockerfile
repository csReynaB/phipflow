FROM rocker/r-ver:4.6.0

ARG QUARTO_VERSION=1.9.37
ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Vienna
ENV RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV MAKEFLAGS="-j4"

# System dependencies for R packages, Quarto, report rendering, fonts, SVG/PDF output,
# GitHub package installation, and compiled packages such as fs, duckdb, stringi, Cairo, svglite.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    build-essential \
    gfortran \
    make \
    cmake \
    pkg-config \
    xz-utils \
    libuv1-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libcairo2-dev \
    libglpk-dev \
    libgmp3-dev \
    libicu-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libgit2-dev \
    libudunits2-dev \
    pandoc \
    fonts-dejavu \
    fonts-montserrat \
    fonts-liberation \
    && rm -rf /var/lib/apt/lists/*

# Install Quarto CLI.
# Quarto bundles its own Pandoc, but system pandoc is useful for knitr/rmarkdown fallback.
RUN wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
      -O /tmp/quarto.deb \
    && apt-get update \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R

WORKDIR /work

# Sanity checks at build time
RUN R --version \
    && quarto --version \
    && pandoc --version \
    && Rscript -e "cat('phiper:', as.character(packageVersion('phiper')), '\n')" \
    && Rscript -e "cat('phiperio:', as.character(packageVersion('phiperio')), '\n')" \
    && Rscript -e "stopifnot(requireNamespace('quarto', quietly = TRUE))"
