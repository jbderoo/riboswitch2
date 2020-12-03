FROM ubuntu

RUN mkdir -p /root/projects/ribo
WORKDIR /root/projects/ribo

RUN apt-get update && \
    apt-get install -y \
    zip unzip \
    tmux \
    python3-venv \
    vim \
    wget

RUN python3 -m venv .venv

RUN /bin/bash -c "source .venv/bin/activate && \
    pip install \
    pandas \
    matplotlib \
    biopython"

RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
    chmod +x Miniconda3-latest-Linux-x86_64.sh  && \
    bash ./Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local  && \
    conda config --add channels defaults  && \
    conda config --add channels bioconda  && \
    conda config --add channels conda-forge  && \
    conda install -q -y --prefix /usr/local locarna

COPY . .

CMD ["/bin/bash"]
