FROM ubuntu:22.04

# Install dependencies including CA certificates and git-lfs
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    build-essential \
    ocl-icd-libopencl1 \
    ocl-icd-opencl-dev \
    clinfo \
    pocl-opencl-icd \
    git \
    git-lfs \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm /tmp/miniconda.sh

ENV PATH=/opt/miniconda/bin:$PATH

RUN conda init bash

# Clone the Git repository without LFS files initially
RUN GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/rajkumar1501/genetic-smiles-optimizer.git /opt/genetic-smiles-optimizer

# Accept conda Terms of Service for required channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create the conda environment from the environment.yml file
RUN conda env create -f /opt/genetic-smiles-optimizer/environment.yml

SHELL ["conda", "run", "-n", "mol2mol_env", "/bin/bash", "-c"]

CMD ["/bin/bash"]
