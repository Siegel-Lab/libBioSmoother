#!/bin/bash

mkdir -p $( dirname $2 )

BUILD_TIME=$( date "+%y-%m-%d-%H:%M:%S" )

# configure the new verison file
NEW_FILE=$(cat $1 | sed "s/@BUILD_TIME@/$BUILD_TIME/g")

echo -n "$NEW_FILE" > $2

# set timestamp to that of original file to hide change from cmake
touch -r $1 $2
