redis-server /usr/local/etc/redis.conf &
REDIS_PID=$!
mongod --config /usr/local/etc/mongod.conf &
rq-dashboard 2> /dev/null &
open "http://0.0.0.0:9181/"

read

killall mongod
kill -9 $REDIS_PID

echo "Terminated"
