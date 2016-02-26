FROM ubuntu:14.04.3
MAINTAINER Tomaz Curk <tomazc@gmail.com>

# thanks to https://github.com/bschiffthaler/ngs/blob/master/base/Dockerfile
# and https://github.com/AveraSD/ngs-docker-star/blob/master/Dockerfile

RUN useradd -m -d /home/icuser icuser

# update system
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    wget \
    g++ \
    zlib1g-dev \
    libzmq-dev \
    make \
    python3 \
    python3-pip \
    python3-setuptools \
    python-virtualenv \
    python-pip

RUN apt-get autoclean -y && \
    apt-get autoremove -y

#################
### RNA-star

# compile STAR from source
WORKDIR /tmp/STAR
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz
RUN tar -xvzf STAR_2.4.2a.tar.gz
WORKDIR /tmp/STAR/STAR-STAR_2.4.2a/source
RUN make STAR
RUN mkdir -p /home/icuser/bin && cp STAR /home/icuser/bin
WORKDIR /tmp
RUN rm -rfv STAR


#################
#### iCount
ADD . /home/icuser/iCount_src

RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser
RUN virtualenv -p python3 /home/icuser/.icountenv

WORKDIR /home/icuser/iCount_src
RUN ../.icountenv/bin/pip install --upgrade -r requirements.txt
RUN ../.icountenv/bin/python setup.py develop

EXPOSE 6543

ENTRYPOINT ["/home/icuser/.icountenv/bin/pserve", "development.ini"]
