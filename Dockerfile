FROM ubuntu:24.04
LABEL maintainer="Yuichi Shiraishi <friend1ws@gmail.com>" 

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt install -y \
    git \
    wget \
    bzip2 \
    make \
    cmake \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    python3 \
    python3-pip

RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
	tar jxvf samtools-1.21.tar.bz2 && \
	cd samtools-1.21 && \
	./configure && make && make install && \
	cd .. && rm -rf samtools-1.21 samtools-1.21.tar.bz2

RUN apt update && apt install -y jellyfish

RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth && \
    chmod +x mosdepth && \
    mv mosdepth /usr/local/bin/

RUN python3 -m pip install --upgrade setuptools --break-system-packages
# RUN python3 -m pip install boto3

RUN git clone https://github.com/friend1ws/ascairn.git && \
	cd ascairn && \
	python3 -m pip install . --break-system-packages && \
	which ascairn 

