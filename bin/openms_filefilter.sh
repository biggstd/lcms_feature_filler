#!/bin/bash

while IFS=, read -r sample_id feature_id rtmin rtmax mzmin mzmax
do
    FileFilter -in $2 -out ${sample_id}__${feature_id}.mzML -mz ${mzmin}:${mzmax} -rt ${rtmin}:${rtmax} 
    collect_features.py --mi-results ${mi_results} --sample-feature ${sample_id}__${feature_id}.mzML
    rm ${sample_id}__${feature_id}.mzML
done < $1