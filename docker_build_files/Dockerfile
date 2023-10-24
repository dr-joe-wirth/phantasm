# syntax=docker/dockerfile:1

# Author: Joseph S. Wirth

# start with ubuntu base
FROM ubuntu:22.04

# install nano, wget, git, python, curl, libcurl, mummer, and blast+
RUN apt-get update
RUN apt-get -y install wget
RUN apt-get -y install python3.9 python3-pip python3-setuptools python3-dev
RUN apt-get -y install curl libcurl4-openssl-dev
RUN apt-get -y install mummer
RUN apt-get -y install ncbi-blast+

# install and update locales (required for R to function properly)
RUN apt-get -y install locales
RUN LC_ALL=C.UTF-8
RUN dpkg-reconfigure locales
RUN echo 'export LC_ALL=C.UTF-8' >> ~/.bashrc

# install R
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y install r-base

# install required R libraries
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gplots',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dendextend',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('DECIPHER')"

# install required python modules
RUN pip3 install --upgrade pip
RUN pip3 install scipy
RUN pip3 install Bio
RUN pip3 install textdistance
RUN pip3 install numpy
RUN pip3 install parasail
RUN pip3 install rpy2
RUN pip3 install pyani

# create the exec directory to hold muscle, fasttree, and iqtree
RUN mkdir /exec

# install fasttree
RUN wget http://www.microbesonline.org/fasttree/FastTreeMP
RUN mv FastTreeMP /exec/FastTreeMP
RUN chmod +x /exec/FastTreeMP

# install muscle
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
RUN mv muscle5.1.linux_intel64 /exec/muscle
RUN chmod +x /exec/muscle

# install iqtree
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz
RUN tar xzf iqtree-1.6.12-Linux.tar.gz
RUN rm iqtree-1.6.12-Linux.tar.gz
RUN mv iqtree-1.6.12-Linux /exec/iqtree

# get xenoGI v3.1.2
RUN wget https://github.com/dr-joe-wirth/xenoGI/archive/refs/tags/v3.1.2.tar.gz
RUN tar xzf v3.1.2.tar.gz
RUN rm v3.1.2.tar.gz

# get phantasm v1.1.2
RUN wget https://github.com/dr-joe-wirth/phantasm/archive/refs/tags/v1.1.2.tar.gz
RUN tar xzf v1.1.2.tar.gz
RUN rm v1.1.2.tar.gz
RUN mv /phantasm-1.1.2 /phantasm

# copy the param file
COPY ./docker_param.py /phantasm/param.py

# modify the variable PHANTASM_PY so that it references the alias
WORKDIR /phantasm
RUN cat phantasm.py | sed -E "s/^(PHANTASM_PY = ).+$/\1'phantasm'/g" > .temp
RUN mv .temp phantasm.py

# modify the phantasm.py executable so that it can be invoked directly from the command line
RUN sed '1s/^/#!\/usr\/bin\/python3\n/' /phantasm/phantasm.py > .temp
RUN mv .temp phantasm
RUN chmod +x /phantasm/phantasm

# add executables to the path
ENV PATH="${PATH}:/phantasm:/exec:/exec/iqtree/bin:"

# create the volume and set it as the working directory
VOLUME ["/data"]
WORKDIR /data

# add an entry point so that phantasm can be called without entering the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
