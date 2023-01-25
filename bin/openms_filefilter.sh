#!/bin/bash

while IFS=, read -r sample_id feature_id rtmin rtmax mzmin mzmax
do
    FileFilter -in $2 -out ${sample_id}__${feature_id}.mzML -mz ${mzmin}:${mzmax} -rt ${rtmin}:${rtmax} 
done < $1