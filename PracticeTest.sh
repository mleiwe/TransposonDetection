#!/bin/bash
forward_read=$1
reverse_read=$2

stem=$(basename "$forward_read" | cut -d_ -f1-2)
echo "$stem"
echo "$stem"_test""