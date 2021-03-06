FROM bmennis/majiq:latest

#Installing Conda 
ENV PATH /opt/conda/lib/R/library/DEXSeq/python_scripts:/opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    conda clean -afy

# Installing Bioconda 
RUN conda config --add channels defaults &&\
conda config --add channels bioconda &&\
conda config --add channels conda-forge 

# Installing Mamba
RUN conda install mamba -n base -c conda-forge -y 

# Installing Samtools
RUN apt update && apt install samtools pigz parallel -y

# Installing Payload
RUN mamba install -c r r-reshape2 -y
RUN mamba install -c conda-forge r-dplyr r-plyr  r-optparse r-feather r-reshape -y
RUN mamba install -c bioconda multiqc htseq bioconductor-dexseq bioconductor-isoformswitchanalyzer bioconductor-rhdf5 bioconductor-deseq2 trimmomatic bioconductor-rsubread kallisto star -y

WORKDIR /MOUNT

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /ENTRYPOINT.sh
COPY [--chown=nobody:nogroup] ./deseq2_dexseq_isoformswitchanalyzer.Rscript /deseq2_dexseq_isoformswitchanalyzer.Rscript
RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]