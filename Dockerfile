FROM ubuntu:22.04

# Install prerequisites + git-lfs configured
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget ca-certificates build-essential ocl-icd-libopencl1 ocl-icd-opencl-dev \
    clinfo pocl-opencl-icd git git-lfs \
  && rm -rf /var/lib/apt/lists/* \
  && git lfs install --system

# Install Miniconda silently
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm /tmp/miniconda.sh && \
    /opt/miniconda/bin/conda clean -afy

ENV PATH=/opt/miniconda/bin:$PATH

# Clone repo and pull LFS
RUN GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/rajkumar1501/genetic-smiles-optimizer.git /opt/genetic-smiles-optimizer && \
    cd /opt/genetic-smiles-optimizer && \
    git lfs pull

# Accept Conda ToS and build environment
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda env create -f /opt/genetic-smiles-optimizer/environment.yml && \
    conda clean -afy

# Ensure entrypoint activates the environment
ENV PATH=/opt/miniconda/envs/mol2mol_env/bin:/opt/miniconda/bin:$PATH

RUN conda init bash && \
    echo "conda activate mol2mol_env" >> ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

CMD ["bash"]
