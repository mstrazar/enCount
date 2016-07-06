# Setup docker on Linux server

- install:
sudo apt-get install nfs-common

- install plug in:
https://github.com/ContainX/docker-volume-netshare

- make it start on boot:
https://github.com/ContainX/docker-volume-netshare/issues/35

- test:
docker run -i -t --volume-driver=nfs -v 192.168.27.50/encode:/endata ubuntu /bin/bash

# Run it in development
docker-compose up

# Run it in production
# (Before running it for the first time, ensure that a data volume for mongo exists:
# run:docker volume create --name mongodata
# )
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d

# check if main sees proper /endata folder:
docker-compose run main /bin/bash

