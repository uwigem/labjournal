#!/usr/bin/env bash

# mkweeks.sh: auto generate weekly journal files

# which week number to start generating on
begin=$(find week-*.md 2> /dev/null | wc -l)
# how many new weeks should be generated
end=${1:-10}

# first day of week 1 in seconds since epoch
firstdate=$(date -d "June 18 2018" +%s)

# loop through weeks for each file
for week in $(seq $begin $end); do
  fname="week-"$(printf "%02d" $week)".md"
  echo -e "# Simulations Week $week\n" > $fname

  # loop through days of week
  for day in {0..6}; do
    curdate=$(($firstdate+86400*((week-1)*7+$day)))
    echo -e "## $(date +'%A, %B %d' -d @$curdate)\n" >> $fname
  done
done
