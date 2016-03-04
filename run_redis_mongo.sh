redis-server /usr/local/etc/redis.conf &
mongod --config /usr/local/etc/mongod.conf &
open "http://0.0.0.0:9181/"
rq-dashboard
