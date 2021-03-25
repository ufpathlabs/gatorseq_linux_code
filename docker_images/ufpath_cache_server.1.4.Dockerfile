#DockerFile for vardict
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.2"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-devel which
ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh
RUN sh Miniconda3-py38_4.8.3-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"
	
RUN conda install -c conda-forge -c anaconda -c bioconda -c r -c kayarre \
	pyyaml==5.1.2 \
	gitpython==3.1.8 \
	openpyxl==3.0.5 \
	seaborn==0.11.0 \
	matplotlib==3.3.2 \
	numpy==1.19.1 \
	pandas==1.1.2 \
	scipy==1.5.2 \
	requests==2.24.0 \
	mysql-connector-python==8.0.18 \
	multiprocess==0.70.10 \
	xmltodict==0.12.0 \
	xmlschema==1.2.5 \
	filelock==3.0.12 \
	awscli==1.18.170
	
	
#RUN conda install -c conda-forge -c anaconda -c bioconda -c r -c kayarre \
#	pyyaml==3.12 \
#	gitpython==2.1.11 \
#	openpyxl==2.4.10 \
#	seaborn==0.11.0 \
#	matplotlib==3.3.2 \
#	numpy==1.19.1 \
#	pandas==1.1.1 \
#	scipy==1.5.2 \
#	requests==2.24.0 \
#	mysql-connector-python==8.0.18 \
#	multiprocess==0.70.10 \
#	xmltodict==0.12.0 \
#	xmlschema==1.2.5
#	
	
	
RUN pip install hl7==0.3.4	
RUN pip install untangle==1.1.0
RUN pip install xlrd
#RUN pip install  rtg-tools
#RUN pip install  hap.py

#RUN rm Miniconda3-4.3.30-Linux-x86_64.sh 


CMD ["bash"]

