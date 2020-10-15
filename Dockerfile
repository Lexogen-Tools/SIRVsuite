FROM python:3.6-buster

RUN pip3 install --upgrade pip \
&& pip3 install numpy>=1.13.3 \
&& pip3 install pycairo>=1.16.0 \
&& pip3 install scipy>=1.0.0 \
&& pip3 install pysam>=0.15.0 \
&& pip3 install pandas>=0.24.0 \
&& pip3 install matplotlib==3.0.0 \
&& pip3 install gtfparse>=1.2.0 \
&& pip3 install pyBigWig>=0.3.17

# TODO: ADD github repo to install the package after switching to public
RUN pip3 install git+http://my_token:sZtLBXmrwFFzmvLiyp-c@10.90.1.56:10080/Bioinfo/sirv-suite.git@SIRVsuite_v0.1
