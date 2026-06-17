# CellWarp -- reproducibility image (Ubuntu 22.04 + Python 3.12).
#
# Ubuntu 22.04 ships Python 3.10; this project pins Python >=3.12,<3.13, so the
# image installs 3.12 from the deadsnakes PPA. The four reproduction gates run
# as a build step, so a successful `docker build` itself certifies that the
# deposit reproduces.
#
#   docker build -t cellwarp .          # build + in-build gate certification
#   docker run  --rm cellwarp           # re-run the four gates + no-download fast-path
#
# This image targets linux/amd64: the deadsnakes PPA -- the validated source of
# Python 3.12 on Ubuntu 22.04 -- ships amd64 only, so the platform is pinned for
# a build that works the same everywhere. On x86-64 reviewer hosts it builds
# natively; on Apple Silicon / ARM it builds under emulation automatically, so
# `docker build` and `docker run` work verbatim on both. The four gates are
# arch-robust: md5 and packet checks are arch-independent and the numeric checks
# carry tolerance.
#
# Scope: the image reproduces the four gates and the no-download fast-path
# (deposited centroids). The full reproduce/run_all.sh pipeline additionally
# needs external CELLxGENE Census data (see DATA_SOURCES.md) and is
# intentionally NOT baked into the image.

FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Python 3.12 via deadsnakes (Ubuntu 22.04's stock python3 is 3.10, which the
# project's requires-python rejects). software-properties-common provides
# add-apt-repository; gnupg + dirmngr are needed to import the PPA signing key
# (they are only *recommended* by software-properties-common, so with
# --no-install-recommends they must be named explicitly, else the key import
# fails); git + ca-certificates round out the toolchain.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        dirmngr \
        git \
        gnupg \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        python3.12 \
        python3.12-venv \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /cellwarp

# Build the venv with Python 3.12; .venv/bin/python is the interpreter every
# gate assumes. Put it on PATH so `python`/`pip` resolve to the venv.
RUN python3.12 -m venv /cellwarp/.venv
ENV PATH="/cellwarp/.venv/bin:${PATH}"

# Bake the deposit source in (the .dockerignore keeps .git and local cruft out).
COPY . /cellwarp

# Install strictly via the existing pyproject pins (wheels; ~1 min, no compiler).
RUN .venv/bin/pip install --upgrade pip \
    && .venv/bin/pip install -e ".[dev]"

# Certify reproduction at build time: a green build == all four gates pass.
RUN echo "### GATE 1: validate.py ###"      && .venv/bin/python reproduce/validate.py \
    && echo "### GATE 2: pytest -q ###"       && .venv/bin/python -m pytest -q \
    && echo "### GATE 3: packet --verify ###" && .venv/bin/python scripts/build_submission_packet.py --verify \
    && echo "### GATE 4: manuscript md5 ###"  && md5sum -c reproduce/MANUSCRIPT_MD5

# Default `docker run` target: re-run the four gates, then the no-download
# fast-path (obs/null ~ 0.522 from deposited centroids, ~2 min, no network).
CMD echo "### GATE 1: validate.py ###"      && .venv/bin/python reproduce/validate.py \
    && echo "### GATE 2: pytest -q ###"       && .venv/bin/python -m pytest -q \
    && echo "### GATE 3: packet --verify ###" && .venv/bin/python scripts/build_submission_packet.py --verify \
    && echo "### GATE 4: manuscript md5 ###"  && md5sum -c reproduce/MANUSCRIPT_MD5 \
    && echo "### fast-path (no download) ###" && .venv/bin/python reproduce/fast_path.py
