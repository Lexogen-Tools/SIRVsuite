FROM python:3.6-buster

# install numpy prior to pyBigWig (otherwise it causes problems)
RUN pip3 install numpy==1.19.2
# install from github
RUN pip3 install git+https://github.com/Lexogen-Tools/SIRVsuite
