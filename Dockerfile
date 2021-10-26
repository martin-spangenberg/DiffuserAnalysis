FROM ubuntu:focal

RUN apt -y update
RUN apt -y install wget tar python3 python3-tk xauth

RUN wget https://root.cern/download/root_v6.24.06.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
RUN tar -xzf root_v6.24.06.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
RUN rm root_v6.24.06.Linux-ubuntu20-x86_64-gcc9.3.tar.gz

RUN apt install -y python3-pip
RUN python3 -m pip install numpy matplotlib pyqt5

COPY scripts/ /scripts/

COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod a+rwx /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

