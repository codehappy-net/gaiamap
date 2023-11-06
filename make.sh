#!/bin/sh
g++ -O3 -I/data/libcodehappy/inc -c gaiamap.cpp -o gaiamap.o
g++ -O3 gaiamap.o /data/libcodehappy/bin/libcodehappy.a -lpthread -o gaiamap
