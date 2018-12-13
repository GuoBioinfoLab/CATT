FROM ubuntu:trusty
MAINTAINER Siyi Chen <chensiyi@hust.edu.cn>

RUN mkdir -p /catt
WORKDIR /catt
COPY tools/ bioTSApypy.py catt.py initialize.py pybam.py  testSample.fq /catt/
COPY resource/ /catt/resource/
RUN apt-get update && apt-get install -y gcc \
                                        g++ \
                                        unzip \
                                        python \
                                        zlib1g-dev \
                                        make \
                                        wget \
   && wget http://144.34.223.36/catt/bowtie2-2.3.4.3-linux-x86_64.zip && wget http://144.34.223.36/catt/samtools-1.9.tar.bz2  && wget http://144.34.223.36/catt/pypy3-v6.0.0-linux64.tar.bz2 \
   && tar -xjf pypy3-v6.0.0-linux64.tar.bz2 && tar -xjf samtools-1.9.tar.bz2 && unzip bowtie2-2.3.4.3-linux-x86_64.zip \
   && cd pypy3-v6.0.0-linux64/bin && cp /catt/get-pip.py ./ && ./pypy3 get-pip.py \
   && ./pip install numpy && ./pip install cython && ./pip install pandas && ./pip install biopython && ./pip install cffi \
   && cd /catt/samtools-1.9 &&  ./configure --without-curses --disable-bz2 --disable-lzma && make && make install \
   && cd /catt && rm pypy3-v6.0.0-linux64.tar.bz2 samtools-1.9.tar.bz2 bowtie2-2.3.4.3-linux-x86_64.zip \
   && ./pypy3-v6.0.0-linux64/bin/pypy3 initialize.py \
   && rm -rf /var/lib/apt/lists/* && apt-get purge -y --auto-remove $buildDeps && rm /var/log/dpkg.log /var/log/alternatives.log /var/log/apt/*.log \
   && chmod +x catt.py

EXPOSE 80