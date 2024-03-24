# Represents docker.io/library/gradle:jdk11-openj9 at 2021-01-14
FROM gradle@sha256:7e07e513b83e6a7790f0cb30820f4142b96ab7ceaac67865965b2127734c2c3d AS java_builder

RUN apt-get update --fix-missing

#---------------------------------------------------------------------------------
# Installing nextflow

ARG NEXTFLOW_REVISION="a96233af6fb505a1716a70408d46218f7d45dc96"

WORKDIR /pgscalculator

RUN wget "https://github.com/nextflow-io/nextflow/archive/${NEXTFLOW_REVISION}.zip" \
  --quiet -O nextflow.zip && \
  unzip -q nextflow.zip && \
  mv "nextflow-${NEXTFLOW_REVISION}" "nextflow" && \
  rm nextflow.zip

WORKDIR /pgscalculator/nextflow

RUN gradle --daemon && BUILD_PACK=1 gradle compile packAll

WORKDIR /pgscalculator


#---------------------------------------------------------------------------------
# Creating runtime image with software built in builder image

#ubuntu:23.04
FROM ubuntu@sha256:04714a1bfbb2d8b5390b5cc0c055e48ebfabd4aa395821b860730ff3277ed74a

COPY --from=java_builder /pgscalculator /pgscalculator

WORKDIR /pgscalculator

RUN mv nextflow/build/releases/nextflow-21.12.1-edge-all /usr/bin/nextflow && \
  chmod +x /usr/bin/nextflow

RUN apt-get update --fix-missing && \
    apt-get install gawk graphviz dos2unix pigz \
    libbz2-dev liblzma-dev vim -y

# Install other dependencies and utilities
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    curl \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


#RUN apt-get update && \
#    apt-get install -y --no-install-recommends ca-certificates && \
#    update-ca-certificates && \
#    wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
#    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
#    rm /tmp/miniconda.sh && \
#    /opt/conda/bin/conda clean -tipsy && \
#    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
#    echo "conda activate base" >> ~/.bashrc
#RUN /opt/conda/bin/conda install -c conda-forge mamba


# Install sbayesR (GCTB version)
WORKDIR /pgscalculator
RUN wget https://cnsgenomics.com/software/gctb/download/gctb_2.05beta_Linux.zip --no-check-certificate && \
  unzip gctb_2.05beta_Linux.zip && \
  mv gctb_2.05beta_Linux/gctb /usr/bin/gctb 

# Install PRS-CS
WORKDIR /repos
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/getian107/PRScs.git PRScs
WORKDIR /pgscalculator

## python 2 micromamba env for PRSCS
# Set the environment variable for micromamba
ENV MAMBA_ROOT_PREFIX=/micromamba
# Ensure the MAMBA_ROOT_PREFIX directory exists
RUN mkdir -p $MAMBA_ROOT_PREFIX
# Copy micromamba from the official image
COPY --from=mambaorg/micromamba:latest /usr/bin/micromamba /usr/local/bin/micromamba
# Initialize micromamba for the nextflow user
RUN micromamba create -n py27 python=2.7 -c conda-forge --yes && \
    micromamba install -n py27 scipy -c conda-forge --yes && \
    micromamba install -n py27 h5py -c conda-forge --yes


## Install python in the right order (ensure python 3.11)
#RUN apt-get update && apt-get install -y software-properties-common
#RUN apt-get update && apt-get install -y python3.11
#RUN apt-get update && apt-get install -y python3.11-dev
#RUN apt-get update && apt-get install -y python3.11-venv
#RUN apt-get update && apt-get install -y python3.11-distutils
#RUN apt-get update && apt-get install -y python3-scipy


# Download and install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/


RUN apt-get update --fix-missing && \
  apt-get install -y openjdk-11-jre

# Download and install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip && \
    unzip plink2_linux_x86_64_20240105.zip && \
    mv plink2 /usr/local/bin/

# Final conf
RUN useradd -rm -s /bin/bash -g users nextflow
#RUN useradd -rm -s /bin/bash -g users -u 1000 nextflow
RUN rm -rf /pgscalculator/* && chown nextflow:users /pgscalculator && chmod g+s /pgscalculator
RUN chmod -R 777 /micromamba

USER nextflow


