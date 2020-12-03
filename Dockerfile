################## BASE IMAGE ######################
FROM biocontainers/biocontainers:latest

################## INSTALLATION ######################

RUN conda install bwa=0.7.17
RUN conda install samtools=0.1.16

WORKDIR /data/

CMD ["bwa"]
CMD ["samtools"]
