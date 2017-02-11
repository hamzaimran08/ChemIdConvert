FROM informaticsmatters/rdkit:Release_2016_09_2
MAINTAINER Daniel Bachler <daniel@douglasconnect.com>

RUN apt-get update && \
    apt-get install -y --no-install-recommends python-pip && \
    rm -rf /var/lib/apt/lists/*
COPY python-flask-server/requirements.txt /python-flask-server/requirements.txt
RUN pip install -r /python-flask-server/requirements.txt --upgrade
COPY python-flask-server/ /python-flask-server/

EXPOSE 8080
CMD ["python", "/python-flask-server/app.py"]
