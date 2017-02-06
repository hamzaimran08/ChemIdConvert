FROM informaticsmatters/rdkit:Release_2016_09_2
MAINTAINER Daniel Bachler <daniel@douglasconnect.com>

RUN apt-get update && \
    apt-get install -y --no-install-recommends python-pip && \
    pip install connexion cirpy && \
    apt-get remove -y python-pip && \
    rm -rf /var/lib/apt/lists/*

COPY python-flask-server/ /python-flask-server/

EXPOSE 8080
CMD ["python", "/python-flask-server/app.py"]
