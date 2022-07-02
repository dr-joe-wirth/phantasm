# syntax=docker/dockerfile:1

# use R installation as the base
FROM r-base:4.1.2

# install nano, git, python, curl, libcurl, and mummer
RUN apt update
RUN apt-get install nano
RUN apt-get -y install git
RUN apt-get install -y python3.9 python3-pip python3-setuptools python3-dev
RUN apt-get install -y curl
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get install -y mummer

# install required R libraries
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gplots',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dendextend',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('DECIPHER')"

# install required python modules
RUN pip3 install scipy
RUN pip3 install Bio
RUN pip3 install textdistance
RUN pip3 install numpy
RUN pip3 install parasail
RUN pip3 install rpy2
RUN pip3 install pyani

# install blast+
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-2.12.0+-x64-linux.tar.gz
RUN mv ncbi-blast-2.12.0+ blast
RUN rm ncbi-blast-2.12.0+-x64-linux.tar.gz

# create the exec directory to hold muscle, fasttree, and iqtree
RUN mkdir /exec

# install fasttree
RUN wget http://www.microbesonline.org/fasttree/FastTreeMP
RUN mv FastTreeMP /exec/FastTreeMP
RUN chmod +x /exec/FastTreeMP

# install muscle
RUN wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
RUN tar xzf muscle3.8.31_i86linux64.tar.gz
RUN rm muscle3.8.31_i86linux64.tar.gz
RUN mv muscle3.8.31_i86linux64 /exec/muscle

# install iqtree
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz
RUN tar xzf iqtree-1.6.12-Linux.tar.gz
RUN rm iqtree-1.6.12-Linux.tar.gz
RUN mv iqtree-1.6.12-Linux /exec/iqtree

# clone xenoGI and phantasm
RUN git clone -b devel https://github.com/ecbush/xenoGI.git
RUN git clone https://github.com/dr-joe-wirth/phantasm.git

# copy the param file
COPY ./docker_param.py /phantasm/param.py

# make an alias for easily calling phantasm
RUN echo "alias phantasm='python3 /phantasm/phantasm.py'" >> ~/.bashrc
