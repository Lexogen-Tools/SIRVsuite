FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

ENV HTSLIB_CONFIGURE_OPTIONS=--enable-plugins

RUN apt-get update && apt-get install -y python3.6
RUN apt-get install -y python3-pip \
&& apt-get install -y libcairo2-dev
RUN pip3 install pandas \
&& pip3 install cython \
&& pip3 install pyranges \
&& pip3 install pycairo \
&& pip3 install biopython \
&& pip3 install pysam

COPY Source /Source
ADD Source /Source/