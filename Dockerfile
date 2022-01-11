FROM ubuntu:18.04

ARG SAM_TOOLS_VERSION="1.3.1"
ARG BCF_TOOLS_VERSION="1.3.1"
ARG HTS_LIB_VERSION="1.3.2"
ENV PYTHONPATH="/rabid_seq"

RUN mkdir /download && \
    mkdir /src && \
    mkdir -p $PYTHONPATH
RUN apt update
RUN apt install -y gcc g++ wget \
    unzip \
    make \
    autoconf \
    automake \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev

RUN wget -O /download/mocha.zip http://software.broadinstitute.org/software/mocha/bio-mocha_$MOCHA_VERSION.zip
RUN wget -O /download/gtc_2_vcf.zip http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_$GTC_2_VCF_VERSION.zip
RUN wget -O /download/samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/$SAM_TOOLS_VERSION/samtools-$SAM_TOOLS_VERSION.tar.bz2
RUN wget -O /download/bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/$BCF_TOOLS_VERSION/bcftools-$BCF_TOOLS_VERSION.tar.bz2
RUN wget -O /download/htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTS_LIB_VERSION/htslib-$HTS_LIB_VERSION.tar.bz2
RUN wget -O $EXECUTION_DIRECTORY/cromwell-$CROMWELL_VERSION.jar https://github.com/broadinstitute/cromwell/releases/download/$CROMWELL_VERSION/cromwell-$CROMWELL_VERSION.jar

RUN cd /download && unzip mocha.zip \
    && mv extendFMT.so /src/extendFMT.so \
    && mv mocha.so /src/mocha.so \
    && mv mocha_plot.R /src/mocha_plot.R \
    && mv mochatools.so /src/mochatools.so \
    && mv pileup_plot.R /src/pileup_plot.R \
    && mv summary_plot.R /src/summary_plot.R \
    && mv trio-phase.so /src/trio-phase.so

RUN cd /download && unzip gtc_2_vcf.zip \
    && mv affy2vcf.so /src/affy2vcf.so \
    && mv gtc2vcf.so /src/gtc2vcf.so \
    && mv gtc2vcf_plot.R /src/gtc2vcf_plot.R

RUN cd /download && tar -xjvf samtools.tar.bz2 \
    && cd samtools-$SAM_TOOLS_VERSION \
    && make \
    && make prefix=/usr/local/bin install \
    && ln -s /usr/local/bin/bin/samtools /usr/bin/samtools

RUN cd /download && tar -xjvf bcftools.tar.bz2 \
    && cd bcftools-$BCF_TOOLS_VERSION \
    && make \
    && make prefix=/usr/local/bin install \
    && ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools

RUN cd /download && tar -xjvf htslib.tar.bz2 \
    && cd htslib-$HTS_LIB_VERSION \
    && make \
    && make install

COPY requirements.txt $PYTHONPATH/requirements.txt
COPY run.sh $PYTHONPATH/mocha/run.sh
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update && apt-get install -y python3-distutils
RUN wget -O /download/get-pip.py https://bootstrap.pypa.io/get-pip.py
RUN cd /download && python3 get-pip.py
RUN cd $PYTHONPATH && pip install -r requirements.txt
