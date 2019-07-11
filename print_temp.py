#!/usr/bin/env python3
import board
import busio
import digitalio
import adafruit_max31865

# -- PT1000 sensor initialization --
# Initialize SPI bus and sensor.
spi = busio.SPI(board.SCK, MOSI=board.MOSI, MISO=board.MISO)
# Chip select of the MAX31865 board.
cs = digitalio.DigitalInOut(board.D5)
# Current setup for pt1000
sensor = adafruit_max31865.MAX31865(spi, cs, wires=3)

print(float(sensor.temperature))
