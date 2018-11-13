FROM ubuntu:16.04

COPY . /fesr

RUN apt-get update
RUN apt-get -y install make
RUN apt-get -y install vim
RUN apt-get -y install gfortran

RUN cd /fesr/fesr_aleph_2015 && make
