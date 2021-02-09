# Select a base image on your server
FROM ubuntu:latest

LABEL maintainer="YJF,2020"

ENV BIO_HOME /usr/local/tool

ENV INSTALL_DIR $BIO_HOME/bin

ENV PATH /usr/local/anaconda3/bin:$PATH

ENV PATH /usr/local/anaconda2/bin:$PATH

ENV PATH $BIO_HOME/circos/bin:$PATH

ENV PATH $BIO_HOME/ViennaRNA/bin:$PATH

ENV PATH $BIO_HOME/RIsearch/bin:$PATH

ENV PATH $BIO_HOME/RNAstructure/exe:$PATH

ENV DATAPATH $BIO_HOME/RNAstructure/data_tables/

ENV PATH $BIO_HOME:$PATH


# Install dependence and NVIDIA-opencl
RUN apt-get update && apt-get install -y vim pciutils wget gcc nvidia-opencl-dev nvidia-opencl-icd-340 


# Copy
COPY uCRISPR $BIO_HOME/uCRISPR


# Wget software
WORKDIR /root

RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh && \
	wget https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh && \
    wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.12.tar.gz && \
    wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureSource.tgz && \
    wget https://rth.dk/resources/risearch/RIsearch-2.1.tar.gz && \
    wget http://circos.ca/distribution/circos-0.69-6.tgz


# Install Anaconda3
RUN bash /root/Anaconda3-2019.07-Linux-x86_64.sh -b -p /usr/local/anaconda3 && \
	rm -f /root/Anaconda3-2019.07-Linux-x86_64.sh
	
RUN /usr/local/anaconda3/bin/pip --no-cache-dir install pyfaidx==0.5.5.2 tensorflow-gpu==1.14.0 dm-sonnet==1.19


# Install Anaconda2
RUN bash /root/Anaconda2-2019.03-Linux-x86_64.sh -b -p /usr/local/anaconda2 && \
	rm -f /root/Anaconda2-2019.03-Linux-x86_64.sh

RUN /usr/local/anaconda2/bin/pip --no-cache-dir install biopython==1.76


# Install Circos
RUN mkdir -p $BIO_HOME/bin && \
	tar -zxvf /root/circos-0.69-6.tgz -C $BIO_HOME && \
	chown -R root:root $BIO_HOME/circos-0.69-6 && \
	ln -s $BIO_HOME/circos-0.69-6 $BIO_HOME/circos && \
	rm -f /root/circos-0.69-6.tgz

RUN mkdir -p /root/.cpan/CPAN/ && \
    apt-get install -y libgd-dev libgd-perl && \
	cpan -i Font::TTF::Font && \
	cpan -i Config::General && \
	cpan -i Clone && \
	cpan -i Math::Bezier && \
	cpan -i List::MoreUtils && \
	cpan -i Regexp::Common && \
	cpan -i Math::Round && \
	cpan -i Math::VecStat && \
	cpan -i Readonly && \
	cpan -i Params::Validate && \
	cpan -i SVG && \
	cpan -i Statistics::Basic && \
	cpan -i Set::IntSpan && \
	cpan -i Text::Format


# Install RIsearch
RUN tar -zxvf /root/RIsearch-2.1.tar.gz -C $BIO_HOME && \
        chown -R root:root $BIO_HOME/RIsearch-2.1 && \
        ln -s $BIO_HOME/RIsearch-2.1 $BIO_HOME/RIsearch && \
        rm -f /root/RIsearch-2.1.tar.gz

# Install ViennaRNA
RUN tar -xzvf /root/ViennaRNA-2.4.12.tar.gz -C $BIO_HOME && \      
        chown -R root:root $BIO_HOME/ViennaRNA-2.4.12 && \
        ln -s $BIO_HOME/ViennaRNA-2.4.12 $BIO_HOME/ViennaRNA && \
        rm -f /root/ViennaRNA-2.4.12.tar.gz

WORKDIR  $BIO_HOME/ViennaRNA

RUN ./configure --prefix=$BIO_HOME/ViennaRNA && make && make install

# Install RNAstructure
RUN tar -xzvf /root/RNAstructureSource.tgz -C $BIO_HOME && \
        chown -R root:root $BIO_HOME/RNAstructure && \
        rm -f /root/RNAstructureSource.tgz 

WORKDIR  $BIO_HOME/RNAstructure 

RUN make all 

# Install uCRISPR
WORKDIR $BIO_HOME/uCRISPR

RUN g++ -o uCRISPR uCRISPR.cpp -std=c++11


EXPOSE 6006

CMD ["/bin/bash"]
