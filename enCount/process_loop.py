import os
import time
import signal

import enCount

# make sure program gets interrupted in a controlled way
stop_it = False


def signal_handler(signal, frame):
    global stop_it
    print('Ctrl+C pressed. Please wait, will stop...')
    stop_it = True
signal.signal(signal.SIGINT, signal_handler)


# call this to empty the redis database
# enCount.queues._redis_conn.flushall()

# remove failed jobs
enCount.queues.failed.empty()

# main scan loop
print('Entering main loop.')

# for testing, download from ENCODE if starting with clean DB
_, online_experiments = enCount.experiments.get_latest()
if not online_experiments:
    print('Checking for new data on ENCODE...')
    online_experiments = enCount.encode.get_online_list()
    enCount.experiments.add_latest_set(online_experiments)

while not stop_it:
    print('Checking for new data on ENCODE...')
    # online_experiments = enCount.encode.get_online_list()

    # compare latest (in DB) and online (on ENCODE) sets of experiments
    # (including their metadata)
    if enCount.experiments.is_updated_set_of_experiments(online_experiments):
        print('An updated list of experiments is available online.')
        print('Will download all new/changed fastq files needed and process '
              'them...')
        print('')
        enCount.experiments.add_latest_set(online_experiments)

    enCount.fastqs.process()
    enCount.experiments.process()
    enCount.mappings.process()
    enCount.queues.print_stats()
#    time.sleep(10)

print('Stopped.')
