#!/usr/bin/env nextflow

// Enable DLS 2.
nextflow.enable.dsl=2


process get_feature_info {
  conda "${params.conda_env}"
  input:
    file sample
    file feature_csv
    file mi_results
  output:
    tuple( path("feature_info__*.csv"), path(sample) )
  script:
  """
  emit_feature_info.py --mi-results '${mi_results}' --sample ${sample} --features '${feature_csv}'
  """
}


process openms_filefilter_feature {
  container = 'docker://docker-registry.phenomenal-h2020.eu/phnmnl/openms'
  input:
    tuple(path(feature_info), path(sample))
    // path "${info[0]}.mzML"
  output:
    file "*.mzML"
  script:
  """
  openms_filefilter.sh '${feature_info}' '${sample}'
  """
}


// while IFS=, read -r sample_id feature_id rtmin rtmax mzmin mzmax
// do
// FileFilter -in !{sample} -out ${sample_id}__${feature_id}.mzML -mz ${mzmin}:${mzmax} -rt ${rtmin]}:${rtmax} 
// done < !{feature_info}


/*******************************************************************************
Extract feature data
--------------------

  Each peak / feature is split into a separate .csv.

*******************************************************************************/
process extract_feature_data {
  conda "${params.conda_env}"
  input:
    each path(sample_feature_file)
    file mi_results
  output:
    stdout emit: summations
    file '*.csv'
  script:
    """
    collect_features.py --mi-results '${mi_results}' --sample-feature '${sample_feature_file}'
    """
}



/*******************************************************************************
Combine peak group CSV files
----------------------------

  The individual .csv files from above are concatenated here. They already
  contain the original filename, so we use nextflow to combine by the peak
  group id names.

*******************************************************************************/
process collate_features {
  // memory '120 GB'
  conda "${params.conda_env}"
  // clusterOptions '--mem=200G'
  publishDir "${params.output_dir}/features", mode: "move", overwrite: true
  input:
    tuple val(feature), file("*.csv")
  output:
    path "${feature}_coll.csv", optional: true
  script:
"""
#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
files = Path('./').glob('*.csv')
if not files:
    exit(0)
dfs = [pd.read_csv(file) for file in files]
df =  pd.concat(dfs, axis=0, ignore_index=True)
df.to_csv("${feature}_coll.csv", index=False)
"""
}




/*******************************************************************************
Plot Peaks
----------

Each selected peak from metaboigniter is plotted and saved.

*******************************************************************************/
// process plot_features {
//   tag "${feature_data.baseName}"
//   conda "${params.conda_env}"
//   publishDir "${params.output_dir}/plots", mode: "move", overwrite: true
//   input:
//     file feature_data
//   output:
//     file "*.png"
//   script:
//     """
//     plot_feature.py --file ${feature_data}
//     """
// }

// // params.conda.enabled  = true
// params.mi_results = "/home/tyer/gits/wsda-smoke-taint/06-feature_analysis/feature_data/AlxVly_01232018.nc"
// params.samples = "${baseDir}/test/samples.csv"
// params.features = "${baseDir}/test/features.csv"
// params.output_dir = "${baseDir}/output"
// params.sample_dir = "/media/tyer/Gamma/wsda_smoke_data/"
// params.conda_env = "/home/tyer/anaconda3/envs/pymc-3.11.4-py39"

workflow {
Channel.fromPath( params.samples )
    .splitText() { it.strip() }
    .map { name -> file("${params.sample_dir}/${name}.mzML") }
    .set { samples }

get_feature_info(samples, file(params.features), file(params.mi_results))
openms_filefilter_feature(get_feature_info.out)
extract_feature_data(openms_filefilter_feature.out.collect().flatten(), file(params.mi_results))
    
extract_feature_data.out[0].collectFile(
  name: 'filled_features.csv', 
  newLine: true, 
  storeDir: "${params.output_dir}")

grouped_features = extract_feature_data.out[1].collect()
  .flatten()
  .map { file -> tuple( file.baseName, file) }
  .groupTuple(by: 0)

collated_features = collate_features(grouped_features)
}
