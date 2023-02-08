#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import argparse
import xarray as xr
import pathlib
import pyopenms
def get_filled_features(feature_id, sample_id, mi_results, sample_exp, 
                        rt_pad_fact=0.0, rt_pad=15, mz_pad_inc=0.000025, ppm_tol=10):
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
    
    intensity_subset = get_subset(
        rt_min, 
        rt_max,
        fs_data['xcmsCamera_mzmin'].values,
        fs_data['xcmsCamera_mzmax'].values, 
        sample_exp)
    
    mz_pad = 0
    while intensity_subset.empty:
        mz_pad += mz_pad_inc
        intensity_subset = get_subset(
            rt_min, 
            rt_max,
            fs_data['xcmsCamera_mzmin'] - mz_pad,
            fs_data['xcmsCamera_mzmax'] + mz_pad, 
            sample_exp)
    
    mz_detla = abs(fs_data['xcmsCamera_mz'].values - intensity_subset['mz'].unique())
    ppm_diff = mz_detla / fs_data['xcmsCamera_mz'].values * 1e6
    
#     print(intensity_subset)
    
    if intensity_subset['mz'].nunique() > 1:
        
        valid_mzs = ppm_diff <= ppm_tol
        
        if not any(valid_mzs):
            print(",".join([sample_id, feature_id, "NaN"]))
            return
        
        
        subset_mzs = intensity_subset['mz'].unique()[valid_mzs]
        subset_mz_diffs = np.abs(subset_mzs - fs_data['xcmsCamera_mz'].values)
        intensity_subset = intensity_subset[intensity_subset['mz'] == subset_mzs[np.argmin(subset_mz_diffs)]]
        
    if intensity_subset.empty:
        print(",".join([sample_id, feature_id, "NaN"]))
        return
        
        
    if all(ppm_diff <= ppm_tol) & (intensity_subset['mz'].nunique() == 1):
        intensity_subset['adjusted_elution_time'] = rt_adj_df['adjusted_elution_time'].iloc[
            new_min_idx:new_min_idx + intensity_subset.shape[0]].values
        intensity_subset['original_elution_time'] = rt_adj_df['original_elution_time'].iloc[
            new_min_idx:new_min_idx + intensity_subset.shape[0]].values
        intensity_subset['Feature'] = feature_id
        intensity_subset['Sample_source'] = sample
        intensity_subset['xcmsCamera_mz'] = fs_data['xcmsCamera_mz'].values
        intensity_subset['xcmsCamera_rt'] = fs_data['xcmsCamera_rt'].values
        intensity_subset['xcmsCamera_rtmin'] = fs_data['xcmsCamera_rtmin'].values
        intensity_subset['xcmsCamera_rtmax'] = fs_data['xcmsCamera_rtmax'].values
        intensity_subset['xcmsCamera_mzmin'] = fs_data['xcmsCamera_mzmin'].values
        intensity_subset['xcmsCamera_mzmax'] = fs_data['xcmsCamera_mzmax'].values
        og_min_rt = rt_adj_df['original_elution_time'][min_idx]
        og_max_rt = rt_adj_df['original_elution_time'][max_idx]
        peak_sum = intensity_subset.loc[intensity_subset['adjusted_elution_time'].between(og_min_rt, og_max_rt), 'inty'].sum()
        print(",".join([sample_id, feature_id, str(peak_sum)]))
        return intensity_subset

    else:
        print(",".join([sample_id, feature_id, "NaN"]))
        return


def get_subset(rtmin, rtmax, mzmin, mzmax, exp):
    data = exp.get2DPeakDataLong(rtmin, rtmax, mzmin, mzmax)
    df = pd.DataFrame(dict(zip(["RT", "mz", "inty"], data)))
    return df


if __name__ == '__main__':
  
    parser = argparse.ArgumentParser()
    parser.add_argument('--mi-results', type=str)
    parser.add_argument('--sample', type=str)
    parser.add_argument('--features', type=str)
    args = parser.parse_args()

    sample_path = str(args.sample)
    sample = str(args.sample).rsplit('.', 1)[0]
    mi_results = xr.open_dataset(args.mi_results)
    features = pd.read_csv(args.features, header=None).values.flatten()

    sample_exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(str(sample_path), sample_exp)

    for feature in features:
        try:
            data = get_filled_features(feature, sample, mi_results, sample_exp)
            if data is not None:
                data.to_csv(f'{feature}.csv', index=False)
        except IndexError as error:
            print(",".join([sample, feature, "IndexError"]))
