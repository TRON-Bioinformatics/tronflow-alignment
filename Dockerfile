FROM nfcore/base:1.9
LABEL authors="Pablo Riesgo Ferreiro" \
      description="Docker image containing all software requirements for the TRON-bioinformatics/tron-bwa pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/tronflow-bwa-1.1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name tronflow-bwa-1.1.0 > tronflow-bwa-1.1.0.yml