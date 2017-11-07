FROM ubuntu:16.04

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    libhighgui-dev \
    libgsl-dev \
    libopencv-core-dev \
    libopencv-highgui-dev \
    libopencv-imgproc-dev \
    libopencv-flann-dev \
    libopencv-photo-dev \
    libopencv-video-dev \
    libopencv-features2d-dev \
    libopencv-objdetect-dev \
    libopencv-calib3d-dev \
    libopencv-ml-dev \
    libopencv-contrib-dev \
    unzip
VOLUME /app
#COPY tracking.zip /tracking.zip
#RUN unzip /tracking.zip
#RUN cd /tracking-paolo-new/CellTracking && make
