FROM jupyter/scipy-notebook
MAINTAINER Daniel Bachler <daniel@douglasconnect.com>

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

# install the RDKit:
RUN conda config --add channels  https://conda.anaconda.org/rdkit
RUN conda install -y nomkl rdkit pandas cairo cairocffi
RUN pip install connexion cirpy

ADD ./python-flask-server/ /python-flask-server/

EXPOSE 8080
ENTRYPOINT ["python", "/python-flask-server/app.py"]


