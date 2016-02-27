import enCount
import time
import signal


stop_it = False

def signal_handler(signal, frame):
    global stop_it
    print('Ctrl+C pressed. Please wait, will stop...')
    stop_it = True

signal.signal(signal.SIGINT, signal_handler)

print('Checking ENCODE for new data')
print('Entering loop')
while not stop_it:
    print("check")
    #enCount.encode.get_list()
    time.sleep(60)

print('Stopped.')
