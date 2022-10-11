#!/bin/bash


GIT_COMMIT_HASH=$(git log -1 --format=%h-%ci | cut -d' ' -f1,2 --output-delimiter -)

GIT_STATUS=$(git status -s)
if [ "${GIT_STATUS}" != "" ]
then
    GIT_STATUS="D-"
fi

CMAKE_VERSION=$2

mkdir -p $( dirname $3 )

# configure the new verison file
NEW_FILE=$(cat $1 | sed "s/@CMAKE_VERSION@/$CMAKE_VERSION/g" \
                  | sed "s/@GIT_COMMIT_HASH@/$GIT_COMMIT_HASH/g" \
                  | sed "s/@GIT_STATUS@/$GIT_STATUS/g")

# sed command does nothing but is required to get the same file format...
OLD_FILE=$(cat $3 | sed "s/@CMAKE_VERSION@/$CMAKE_VERSION/g")

# only write the file if the version actually changed
if [ "$NEW_FILE" != "$OLD_FILE" ]
then
    echo "Writing new version file"
    echo -n "$NEW_FILE" > $3
fi