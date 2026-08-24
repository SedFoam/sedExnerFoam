FROM cbonamy/sedfoam_2512_ubuntu:develop
LABEL maintainer="cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10
ARG DEBIAN_FRONTEND=noninteractive

USER root

WORKDIR /root/
RUN git clone --branch develop --recurse-submodules https://github.com/sedfoam/sedexnerfoam.git --depth 1

WORKDIR /root/sedexnerfoam
RUN /bin/bash -c 'shopt -s expand_aliases && source /openfoam/bash.rc && export FOAM_USER_LIBBIN=$FOAM_SITE_LIBBIN && ./Allwmake -prefix=group'

USER 1001

WORKDIR /home/sudofoam
RUN sudo cp -r /root/sedexnerfoam sedexnerfoam && sudo chown -R sudofoam: /home/sudofoam && sudo rm -rf /root/*

WORKDIR /home/sudofoam
ENV HOME=/home/sudofoam
ENTRYPOINT ["/bin/bash", "-c"]
