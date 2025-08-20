# Represents docker.io/library/gradle:jdk11-openj9 at 2021-01-14
FROM gradle:8.5-jdk11 AS java_builder

RUN apt-get update --fix-missing

#---------------------------------------------------------------------------------
# Installing nextflow using pre-built binary

ARG NEXTFLOW_VERSION="24.10.4"

WORKDIR /pgscalculator

# Download and install Nextflow standalone distribution
RUN curl -fsSL https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-dist > /usr/local/bin/nextflow && \
    chmod +x /usr/local/bin/nextflow && \
    # Pre-download basic plugin dependencies
    nextflow plugin install nf-validation

# Create nextflow cache directory
RUN mkdir -p /root/.nextflow

# Ensure the cache is preserved in the final image
VOLUME /root/.nextflow

#---------------------------------------------------------------------------------
# Compiling b3sum since it is not packaged for debian
FROM rust:1.84-slim-bookworm AS rust_builder

RUN cargo install --version 0.3.8 b3sum

#---------------------------------------------------------------------------------
# Creating runtime image with software built in builder image
FROM eclipse-temurin:11

COPY --from=java_builder /usr/local/bin/nextflow /usr/local/bin/nextflow
COPY --from=java_builder /root/.nextflow /root/.nextflow
COPY --from=rust_builder /usr/local/cargo/bin/b3sum /usr/bin/b3sum

WORKDIR /pgscalculator

# Rest of your installations
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    curl \
    gawk \
    graphviz \
    dos2unix \
    pigz \
    libbz2-dev \
    liblzma-dev \
    vim \
    parallel \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgomp1 \
    git \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('data.table', repos='https://cloud.r-project.org/')"

# Install PRS-CS
WORKDIR /repos
RUN git clone https://github.com/getian107/PRScs.git PRScs

# Install GCTB (sbayesR)
WORKDIR /pgscalculator
RUN wget https://cnsgenomics.com/software/gctb/download/gctb_2.05beta_Linux.zip --no-check-certificate && \
    unzip gctb_2.05beta_Linux.zip && \
    mv gctb_2.05beta_Linux/gctb /usr/bin/gctb 

# Install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/

# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip && \
    unzip plink2_linux_x86_64_20240105.zip && \
    mv plink2 /usr/local/bin/

# Set up directories with correct permissions
RUN useradd -rm -s /bin/bash -g users nextflow && \
    mkdir -p /pgscalculator/.nextflow && \
    chown -R nextflow:users /pgscalculator && \
    chmod -R 777 /pgscalculator && \
    chmod g+s /pgscalculator

# Then set nextflow for offline mode
ENV NXF_OFFLINE=true

USER nextflow

## Initialize Nextflow in the final image
#RUN nextflow -version && \
#    nextflow info && \
#    # Run a simple workflow to ensure dependencies are downloaded
#    echo "workflow { println 'Hello' }" > test.nf && \
#    nextflow run test.nf && \
#    rm test.nf


