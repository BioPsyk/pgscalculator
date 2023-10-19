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

# Represents docker.io/library/adoptopenjdk:11-jre-openj9 at 2021-01-14
FROM adoptopenjdk@sha256:477d0c53aca999692d2432e529af1f7abd715205fcfc36534ac9ff490f4da0e8

COPY --from=java_builder /pgscalculator /pgscalculator

WORKDIR /pgscalculator

RUN mv nextflow/build/releases/nextflow-21.12.1-edge-all /usr/bin/nextflow && \
  chmod +x /usr/bin/nextflow

RUN apt-get update --fix-missing && \
    apt-get install python3 python3-pip python gawk graphviz dos2unix pigz \
    libbz2-dev liblzma-dev vim -y

# Install other dependencies and utilities
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip && \
    unzip plink_linux_x86_64_20231018.zip && \
    mv plink /usr/local/bin/

# Install sbayesR, Install R and dependencies
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('sbayesR', repos='https://cloud.r-project.org/')"

# Install PRS-CS
RUN wget https://github.com/getian107/PRScs/raw/master/PRScs.py -O /usr/local/bin/PRScs.py && \
    chmod +x /usr/local/bin/PRScs.py


# Final conf
RUN useradd -rm -s /bin/bash -g users -u 1000 nextflow
RUN rm -rf /pgscalculator/* && chown nextflow:users /pgscalculator && chmod g+s /pgscalculator

USER nextflow

