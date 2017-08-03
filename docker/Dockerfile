# Docker file that installs docker container for Selenzy
#
# rename this file to "Dockerfile"
# build with: "sudo docker build -t selenzy ."

FROM continuumio/anaconda3 

# Install rdkit
RUN conda install -c rdkit rdkit
RUN conda install -c conda-forge flask-restful
RUN conda install -c anaconda biopython
RUN conda install -c bioconda emboss
RUN conda install -c biobuilds t-coffee

# To be replace by a git clone
RUN wget http://130.88.113.226/selenzy/selenzy.tar.gz
RUN tar -xzvf selenzy.tar.gz

ENTRYPOINT ["python"]

CMD ["/selenzyPro/flaskform.py", "/selenzyPro/uploads", "/selenzyPro/data"]

EXPOSE 5000

# sudo docker run --name selenzy -d 5555:5000 -e LD_LIBRARY_PATH='/opt/conda/bin/../lib' -it /bin/bash
