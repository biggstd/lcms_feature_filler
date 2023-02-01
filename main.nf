#!/usr/bin/env nextflow

// Enable DLS 2.
nextflow.enable.dsl=2


/*******************************************************************************
Extract feature data
--------------------

  Each peak / feature is split into a separate .csv.

*******************************************************************************/
process extract_feature_data {
  conda "${params.conda_env}"
  input:
    file sample_file
    file features
    file mi_results
  output:
    stdout emit: summations
    file '*.csv'
  script:
    """
    collect_features.py --mi-results '${mi_results}' --sample '${sample_file}' \
    --features '${features}'
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


workflow {
Channel.fromPath( params.samples )
    .splitText() { it.strip() }
    .map { name -> file("${params.sample_dir}/${name}.mzML") }
    .set { samples }

  extract_feature_data(samples, file ( params.features ), file ( params.mi_results ) )

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
