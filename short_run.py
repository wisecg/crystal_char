#!/usr/bin/env python3
import time
import datetime as datetime
import subprocess as sp
import board
import busio
import digitalio
import adafruit_max31865

def main():
    # -- PT1000 sensor initialization --
    # Initialize SPI bus and sensor.
    spi = busio.SPI(board.SCK, MOSI=board.MOSI, MISO=board.MISO)
    # Chip select of the MAX31865 board.
    cs = digitalio.DigitalInOut(board.D5)
    # Current setup for pt1000
    sensor = adafruit_max31865.MAX31865(spi, cs, wires=3)

    # -- set parameters --
    #runinfo = open("runinfo.txt", "r")
    #sn = runinfo.read()
    #runinfo.close()

    run_number = input("Run number: ")
    filename = 'run%s_temperature_data.txt' % (run_number)
    start_time = time.time()
    print('Start time: %s' % (time.ctime(start_time)))
    print('Estimated end time: %s' % (time.ctime(start_time + 600)))
    save_path = '/home/pi/tempstudy/data/'
    file = save_path + filename
    f = open(file,'w')
    f.close()
    data = get_temperatures(sensor, file)
    f.close()
    sync_data(filename)

def get_temperatures(sensor, file):
    """
    10 minutes
    600 seconds
    """
    duration = 660
    interval = 1
    temperature_data = []

    for pos in range(1,duration):
        curr_temp = float(sensor.temperature)
        curr_time = time.time()
        print('---------------------------------------------')
        print('Time: %s' % (time.ctime(curr_time)))
        print('Temperature: %f' % (curr_temp))
        f = open(file, 'a')
        f.write('%f %f\n' % (curr_time, curr_temp))
        f.close()
        curr = [curr_time, curr_temp]
        temperature_data.append(curr)
        time.sleep(interval)

    return temperature_data

def write_temperature(file, data):
    f = open(file,'w')
    for i in range(len(data) - 1):
        str = "%f %f\n" % (data[i][0], data[i][1])
        f.write(str)
    f.close()

def sync_data(filename):
    destination = 'ccenpa@d-172-25-100-73.dhcp4.washington.edu'
    location = '~/Desktop/coherent/Analysis/temp_study/temperature_data'
    sp.Popen('scp %s %s:%s' % (filename, destination, location), shell=True)
    print('%s synced to %s' % (filename, location))

if __name__=="__main__":
    main()
