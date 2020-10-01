FROM python:3.6-buster

# install numpy prior to pyBigWig (otherwise it causes problems)
RUN pip3 install numpy==1.19.2
# install from gitlab
RUN pip3 install git+http://my_token:sZtLBXmrwFFzmvLiyp-c@10.90.1.56:10080/Bioinfo/sirv-suite.git@SIRVsuite_v0.1
