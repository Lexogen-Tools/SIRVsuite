FROM python:3.6-buster

# install from gitlab
RUN pip3 install git+http://my_token:sZtLBXmrwFFzmvLiyp-c@10.90.1.56:10080/Bioinfo/sirv-suite.git@SIRVsuite_v0.1
