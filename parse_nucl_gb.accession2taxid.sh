#! /usr/bin/env bash

p="{"
while read -r a v t g
do
   echo -n $p
   p="\"$v\": $t, "
done < <(zcat "${1:-/dev/stdin}" | tail -n +2)
p=${p%,*}
echo -n "$p}"
