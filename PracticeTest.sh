#!/bin/bash
echo "This is a test"
ChkFile=$1
echo $ChkFile
if [ -f "$ChkFile" ]; then
    echo "It works"
else
    echo "It didn't find the file"
fi