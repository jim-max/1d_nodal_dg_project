FROM debian:stable-slim

RUN apt-get update && apt-get install -y \
    clang-tidy \
    clang-format \
    clangd \
    git \
    g++ \
    cmake \
    pre-commit \
    doxygen \
    gnuplot \
    libomp-dev \
    libboost-program-options-dev \
    libboost-thread-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
