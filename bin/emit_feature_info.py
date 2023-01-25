#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import argparse
import xarray as xr
import pathlib
import pyopenms



def emit(mi_results, sample_id, feature_id, mz_pad=0.05, rt_pad_fact=0.0, rt_pad=15):
    rt_adj_df = pd.DataFrame(
    {'adjusted_elution_time': mi_results.sel(
        Sample=sample_id)['adjusted_elution_time'].values,
     'original_elution_time': mi_results.sel(
        Sample=sample_id)['original_elution_time'].values})
    fs_data = mi_results.sel(Sample=sample_id, Feature=feature_id)

    min_idx = (rt_adj_df['adjusted_elution_time'] 
               - fs_data['xcmsCamera_rtmin'].values).abs().idxmin()
    
    max_idx = (rt_adj_df['adjusted_elution_time'] 
               - fs_data['xcmsCamera_rtmax'].values).abs().idxmin()
    
    idx_width = max_idx - min_idx
    new_min_idx = int(min_idx - idx_width * rt_pad_fact - rt_pad) + 1
    new_max_idx = int(max_idx + idx_width * rt_pad_fact + rt_pad)
    
    idx_size = rt_adj_df['adjusted_elution_time'].values.shape[0]
    
    if new_min_idx < 0:
        new_min_idx = 0
        
    if new_max_idx >= idx_size:
        new_max_idx = idx_size - 1
    
    rt_min = rt_adj_df['original_elution_time'].values[new_min_idx]
    rt_max = rt_adj_df['original_elution_time'].values[new_max_idx]
    
    return [
        sample_id,
        feature_id,
        rt_min,
        rt_max,
        fs_data['xcmsCamera_mzmin'].values - mz_pad,
        fs_data['xcmsCamera_mzmax'].values + mz_pad,
    ]
    # print(sample_id)
    # print(feature_id)
    # print(rt_min)
    # print(rt_max)
    # print(fs_data['xcmsCamera_mzmin'].values - mz_pad)
    # print(fs_data['xcmsCamera_mzmax'].values + mz_pad)
          

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mi-results', type=str)
    parser.add_argument('--sample', type=str)
    parser.add_argument('--features', type=str)
    args = parser.parse_args()

    mi_results = xr.open_dataset(args.mi_results)
    sample_path = str(args.sample)
    sample = pathlib.Path(str(args.sample)).name.rsplit('.', 1)[0]
    features = pd.read_csv(args.features, header=None).values.flatten()

    info = pd.DataFrame([emit(mi_results, sample, feature) for feature in features])
    info.to_csv(f"feature_info__{sample}.csv", index=None, header=None)
    
    # for feature in features:
        # emit(mi_results, sample, feature)