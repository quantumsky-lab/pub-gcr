# Debian GNU/Linux 10
FROM python:3 
USER root
RUN mkdir -p /root/data
WORKDIR /root

RUN pip install --ignore-installed --no-cache-dir \
    setuptools \
    ipyparallel \
    argparse \
    scipy \
    google-cloud-storage \
    biopython

RUN pip install pandas

ENV SHELL=/bin/bash

# get JAVA for trimmomatic
RUN apt update && apt install -y default-jdk

# get blast
RUN apt-get update && apt-get install -y ncbi-blast+

# get trimmomatic
COPY ./trimmomatic-0.39.jar /root/

# copy scripts
COPY ./mpileup2matrix.py /root/

# run mpileup2matrix.py
ENTRYPOINT ["python", "-W ignore", "/root/mpileup2matrix.py"]