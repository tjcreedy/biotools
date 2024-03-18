#! /usr/bin/env bash

echo -n "{"

while read -r a v t g
do
   echo -n "\"$v\": $t,"
done < <(zcat "${1:-/dev/stdin}" | tail -n +2)

echo -ne "\b}"
