FROM ubuntu:18.04
RUN ["/bin/bash", "-c", "apt-get update -y --no-install-recommends"]
RUN ["/bin/bash", "-c", "apt-get upgrade -y --no-install-recommends"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends bc"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends r-base"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends samtools"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends bowtie2"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends bedtools"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends parallel"]
RUN ["/bin/bash", "-c", "apt-get install -y --no-install-recommends libgsl-dev"]
RUN ["/bin/bash", "-c", "ln /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgsl.so.0"]
RUN mkdir /home/src
RUN mkdir /home/genomes
RUN mkdir /home/input
RUN mkdir /home/ouput
COPY Genomes/ /home/genomes
COPY HiCUP-master/ /home/src/HiCUP-master
COPY preseq_v2.0/preseq /home/src/preseq
COPY main.sh /home/src/
COPY script_raw.sh /home/src/script_raw.sh
COPY script_bam.sh /home/src/script_bam.sh
COPY script_bam_fast.sh /home/src/script_bam_fast.sh
COPY scriptsR/script_20th_perc.R /home/src/
COPY scriptsR/script_equa.R /home/src
COPY scriptsR/script_predict.R /home/src
COPY scriptsR/script_predict_no_preseq.R /home/src
COPY scriptsR/script_predict_bam.R /home/src
COPY templates/ /home/templates
ENTRYPOINT ["/home/src/main.sh"]
CMD []
