import tensorflow
from tensorflow.python.client import device_lib

for device in device_lib.list_local_devices():
    print(device.name)


