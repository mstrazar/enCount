## Setup docker on Linux server

1. Install NFS:

        sudo apt-get install nfs-common

2. Install plugin for NFS shared volumes, download it from [github](https://github.com/ContainX/docker-volume-netshare).

3. Make the plugin start on boot, follow [instructions](https://github.com/ContainX/docker-volume-netshare/issues/35).

4. Test and check if NFS volume is available:

        docker run -i -t --volume-driver=nfs -v 192.168.27.50/encode:/endata ubuntu /bin/bash

## Run it in development
Create two volumes needed for testing:

       docker volume create --name=mongodata_test
       docker volume create --name=endata_test

Go into the enCount folder and start:

    docker-compose up

## Run it in production
Before running it for the first time, ensure that a data volume for mongo exists:

    run:docker volume create --name mongodata

Then run it by applying the production specification of volumes.

    docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d

Check if main sees proper /endata folder:

    docker-compose run main /bin/bash
