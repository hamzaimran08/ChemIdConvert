Chemical ID Conversion Tool
==========================

This is a service that provides a REST Api to convert
between different chemical notations / Ids.

Running the service
-------------------

The easiest way to run the service is to use docker. Either
build the container yourself like this:

##### Setup the docker container

```
docker build -t chemidconvert .
```

Or get the latest version from our douglas connect repository
(if you have configured gcloud and have access to our docker
registry):

```
gcloud docker pull eu.gcr.io/douglasconnect-docker/chemidconvert:latest
```

##### Running the docker container

```
docker run -it -p 8080:8080 eu.gcr.io/douglasconnect-docker/chemidconvert
```

Using the Api
-------------

To see what operations are currently supported, consult the swagger api file.
The url will then be: (HOSTNAME is either localhost if running docker locally
on linux or the ip of your docker-machine if running locally on windows/osx)

```
http://HOSTNAME:8080/v1/smiles/toInchi?smiles=CCO
```

