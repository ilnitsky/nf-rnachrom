from functools import partial
from math import isnan
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.interpolate as si
import scipy.stats as ss
from scipy.stats._stats_mstats_common import LinregressResult
from tqdm.contrib.concurrent import process_map

from ..api.binops import calculate_dists_to_centers, calculate_dists_to_points
from ..api.formats import Rdc
from ..api.mp import adjust_chunksize
from ..api.schemas import GeneCoord, SplineResult


NO_LINREG_RESULT = LinregressResult(slope=np.nan,
                                    intercept=np.nan,
                                    rvalue=np.nan,
                                    pvalue=np.nan,
                                    stderr=np.nan,
                                    intercept_stderr=np.nan)
NO_SPLINE_RESULT = SplineResult(t=np.array([]), c=np.array([]), k=np.nan)


def _get_cis_coverage_single(rdc_data: Rdc, rna_name: str, gene_coord: GeneCoord) -> pd.DataFrame:
    cis_coverage = rdc_data.read_pixels(rna_name,
                                        value_fields=['signal_count',
                                                      'bg_count',
                                                      'signal_prob',
                                                      'bg_prob'],
                                        chrom_type='cis')

    undefined = (cis_coverage['signal_count'] == 0) | (cis_coverage['bg_count'] == 0)
    cis_coverage['fc'] = cis_coverage['signal_prob'] / cis_coverage['bg_prob']
    cis_coverage = cis_coverage.loc[~undefined, :].reset_index(drop=True)

    bin_centers = (cis_coverage['start'] + cis_coverage['end']) // 2
    dists_from_gene = calculate_dists_to_points(bin_centers, gene_coord.start, gene_coord.end)

    cis_coverage['log_gene_dist'] = np.log10(dists_from_gene)
    cis_coverage['log_fc'] = np.log10(cis_coverage['fc'])

    cis_coverage = cis_coverage.query('signal_count > 0 and bg_count > 0')[['log_gene_dist', 'log_fc', 'signal_count', 'bg_prob']].reset_index(drop=True)

    cis_coverage['chrom'] = gene_coord.chrom
    cis_coverage['rna'] = rna_name
    return cis_coverage


def get_cis_coverage(rdc_data: Rdc, n_cores: int = 1, chunksize: int = 50) -> pd.DataFrame:
    annotation = rdc_data.annotation
    rna_names = list(annotation.keys())
    gene_coords = [annotation[rna_name] for rna_name in rna_names]
    all_chroms = list({item.chrom for item in gene_coords})

    func = partial(_get_cis_coverage_single, rdc_data)
    chunksize = adjust_chunksize(len(rna_names), n_cores, chunksize)
    cis_results = process_map(func,
                              rna_names,
                              gene_coords,
                              max_workers=n_cores,
                              chunksize=chunksize,
                              desc='Getting cis coverage',
                              unit='RNA')

    cis_coverage = pd.concat(cis_results, ignore_index=True)
    cis_coverage['chrom'] = cis_coverage['chrom'].astype(pd.CategoricalDtype(all_chroms))
    cis_coverage['rna'] = cis_coverage['rna'].astype(pd.CategoricalDtype(rna_names))
    return cis_coverage


def _get_chrom_scaling_single(chrom_name: str,
                              chrom_data: pd.DataFrame,
                              degree: int = 3) -> Tuple[str, Tuple[LinregressResult, SplineResult]]:
    chrom_data = chrom_data.sort_values('log_gene_dist', ignore_index=True)
    x = chrom_data['log_gene_dist'].values
    y = chrom_data['log_fc'].values
    try:
        linreg_results: LinregressResult = ss.linregress(x, y)
    except ValueError:
        linreg_results = NO_LINREG_RESULT

    if len(chrom_data) > degree:
        spl = si.splrep(x,
                        y,
                        w=np.ones(x.shape),
                        k=degree)
        spline_results = SplineResult(t=spl[0], c=spl[1], k=spl[2])
    else:
        spline_results = NO_SPLINE_RESULT
    return chrom_name, (linreg_results, spline_results)


def _get_chrom_scaling(cis_coverage: pd.DataFrame,
                       only_fittable: bool = False,
                       degree: int = 3,
                       n_cores: int = 1,
                       chunksize: int = 1) -> Tuple[Dict[str, LinregressResult], Dict[str, SplineResult]]:
    if only_fittable:
        chrom_counts = cis_coverage['chrom'].value_counts()
        fittable_chroms = chrom_counts[chrom_counts > degree].index
        non_fittable_chroms = chrom_counts[chrom_counts <= degree].index
        cis_coverage = cis_coverage.query('chrom in @fittable_chroms').reset_index(drop=True)

    func = partial(_get_chrom_scaling_single, degree=degree)
    chunksize = adjust_chunksize(len(fittable_chroms), n_cores, chunksize)

    results = process_map(func,
                          *zip(*cis_coverage.groupby('chrom')),
                          max_workers=n_cores,
                          desc='Calculating chromosome-level scaling',
                          unit='chrom',
                          chunksize=chunksize)

    linreg_results = {item[0]: item[1][0] for item in results}
    linreg_results.update({chrom: NO_LINREG_RESULT
                           for chrom in non_fittable_chroms})

    spline_results = {item[0]: item[1][1] for item in results}
    spline_results.update({chrom: NO_SPLINE_RESULT
                           for chrom in non_fittable_chroms})

    return linreg_results, spline_results


def _get_rna_scaling_single(rna_name: str,
                            rna_data: pd.DataFrame,
                            degree: int = 3,
                            no_refine: bool = False,
                            max_threshold: float = 0.05) -> Tuple[str, Optional[SplineResult]]:
    rna_data = rna_data.sort_values('log_gene_dist', ignore_index=True)
    x = rna_data['log_gene_dist'].values
    y = rna_data['log_fc'].values
    if len(x) <= degree:
        return rna_name, None
    else:
        spl = si.splrep(x, y, k=degree, w=np.ones(x.shape))
        if any(isnan(item) for record in spl[:2] for item in record):
            return rna_name, None

        if not no_refine:
            scaling_factors = 10 ** si.splev(x, spl)
            cis_scaled_prob = rna_data['bg_prob'] * scaling_factors
            rna_data['scaled_bg_prob'] = cis_scaled_prob / cis_scaled_prob.sum()
            n_total = rna_data['signal_count'].sum()
            pvals = rna_data.apply(lambda row: ss.binomtest(k=row['signal_count'],
                                                            n=n_total,
                                                            p=row['scaled_bg_prob'],
                                                            alternative='greater').pvalue,
                                   axis=1)
            threshold = min(1 / len(x), max_threshold)
            rna_data['outlier'] = (pvals < threshold)
            cleaned_rna_data = rna_data.query('~outlier')

            if cleaned_rna_data.shape[0] > degree:
                x = cleaned_rna_data['log_gene_dist'].values
                y = cleaned_rna_data['log_fc'].values
                spl = si.splrep(x, y, k=degree, w=np.ones(x.shape))

        spline_results = SplineResult(t=spl[0], c=spl[1], k=spl[2])
        return rna_name, spline_results


def _get_rna_scaling(cis_coverage: pd.DataFrame,
                     degree: int = 3,
                     no_refine: bool = False,
                     max_threshold: float = 0.05,
                     n_cores: int = 1,
                     chunksize: int = 50) -> pd.DataFrame:
    func = partial(_get_rna_scaling_single, degree=degree, no_refine=no_refine, max_threshold=max_threshold)
    chunksize = adjust_chunksize(cis_coverage['rna'].nunique(), n_cores, chunksize)
    results = process_map(func,
                          *zip(*cis_coverage.groupby('rna')),
                          max_workers=n_cores,
                          chunksize=chunksize,
                          desc='Calculating RNA-level scaling',
                          unit='RNA')
    return dict(results)


def _refine_rna_splines(rna_splines: Dict[str, SplineResult],
                        chrom_splines: Dict[str, SplineResult],
                        annotation: Dict[str, GeneCoord]) -> Dict[str, SplineResult]:
    rna_splines_refined = dict()
    for rna_name, rna_spline in rna_splines.items():
        if rna_spline is None:
            rna_spline = chrom_splines.get(annotation[rna_name].chrom,
                                           SplineResult(t=np.array([]), c=np.array([]), k=np.nan))
        rna_splines_refined[rna_name] = rna_spline
    return rna_splines_refined


def _rescale_rdc_data_single(rna_name: str,
                             gene_coord: GeneCoord,
                             n_cis_contacts: int,
                             n_trans_contacts: int,
                             rdc_data: Rdc,
                             fill_value: Union[int, float]) -> Tuple[str, pd.DataFrame]:
    cis_pixels = rdc_data.read_pixels(rna_name, value_fields=['raw_bg_prob', 'bg_prob', 'signal_prob'], chrom_type='cis')
    trans_pixels = rdc_data.read_pixels(rna_name, value_fields=['raw_bg_prob', 'bg_prob', 'signal_prob'], chrom_type='trans')
    single_spline = tuple(rdc_data.read_scaling(rna_name))
    gene_start = gene_coord.start
    gene_end = gene_coord.end
    if not isnan(single_spline[2]):
        spl = single_spline
        rel_cis_bin_centers = calculate_dists_to_centers(cis_pixels, gene_start, gene_end)
        rel_cis_bin_starts = calculate_dists_to_points(cis_pixels['start'], gene_start, gene_end)
        rel_cis_bin_ends = calculate_dists_to_points(cis_pixels['end'], gene_start, gene_end)
        scaling_factors = 10 ** ((si.splev(np.log10(rel_cis_bin_starts + 1), spl, ext=3)
                                 + si.splev(np.log10(rel_cis_bin_ends + 1), spl, ext=3)
                                 + si.splev(np.log10(rel_cis_bin_centers + 1), spl, ext=3)) / 3)
        cis_pixels['scaling_factor'] = scaling_factors
    else:
        cis_pixels['scaling_factor'] = 1

    trans_pixels['scaling_factor'] = 1

    cis_scaled_bg_probs = cis_pixels['raw_bg_prob'] * cis_pixels['scaling_factor']
    cis_pixels['bg_prob'] = cis_scaled_bg_probs / cis_scaled_bg_probs.sum()

    trans_scaled_bg_probs = trans_pixels['raw_bg_prob'] * trans_pixels['scaling_factor']
    trans_pixels['bg_prob'] = trans_scaled_bg_probs / trans_scaled_bg_probs.sum()
    total_interested_contacts = n_cis_contacts + n_trans_contacts
    cis_share = n_cis_contacts / total_interested_contacts
    trans_share = n_trans_contacts / total_interested_contacts

    trans_pixels['raw_bg_prob'] *= trans_share
    trans_pixels['bg_prob'] *= trans_share
    trans_pixels['signal_prob'] *= trans_share

    cis_pixels['raw_bg_prob'] *= cis_share
    cis_pixels['bg_prob'] *= cis_share
    cis_pixels['signal_prob'] *= cis_share

    trans_pixels['fc'] = (trans_pixels['signal_prob'] / trans_pixels['bg_prob']).fillna(fill_value)
    cis_pixels['fc'] = (cis_pixels['signal_prob'] / cis_pixels['bg_prob']).fillna(fill_value)
    changed_cols = ['raw_bg_prob', 'bg_prob', 'signal_prob', 'scaling_factor', 'fc']
    all_pixels = pd.concat([trans_pixels, cis_pixels], ignore_index=True)[['chrom'] + changed_cols]
    return rna_name, all_pixels


def _rescale_rdc_data(rdc_data: Rdc,
                      fill_value: Union[int, float] = 1,
                      n_cores: int = 1,
                      chunksize: int = 50):
    cis_contacts_nums = rdc_data.read_rna_attribute_batch('cis_contacts')
    trans_contacts_nums = rdc_data.read_rna_attribute_batch('trans_contacts')
    annotation = rdc_data.annotation
    if not rdc_data.is_scaling_fitted:
        raise Exception

    func = partial(_rescale_rdc_data_single, rdc_data=rdc_data, fill_value=fill_value)
    rna_names = list(annotation.keys())
    gene_coords = (annotation[rna_name] for rna_name in rna_names)
    n_cis_contacts_gen = (cis_contacts_nums[rna_name] for rna_name in rna_names)
    n_trans_contacts_gen = (trans_contacts_nums[rna_name] for rna_name in rna_names)

    chunksize = adjust_chunksize(len(rna_names), n_cores, chunksize)
    results = dict(process_map(func,
                               rna_names,
                               gene_coords,
                               n_cis_contacts_gen,
                               n_trans_contacts_gen,
                               max_workers=n_cores,
                               chunksize=chunksize,
                               desc='Rescaling rdc pixels',
                               unit='RNA'))
    changed_cols = ['raw_bg_prob', 'bg_prob', 'signal_prob', 'scaling_factor', 'fc']
    for col in changed_cols:
        rdc_data.write_pixels_column_batch(col, results)


def calculate_scaling_splines(rdc_data: Rdc,
                              degree: int = 3,
                              no_refine: bool = False,
                              max_threshold: float = 0.05,
                              fill_value: Union[int, float] = 1,
                              n_cores: int = 1) -> None:
    cis_coverage = get_cis_coverage(rdc_data)
    _, chrom_spline_scaling = _get_chrom_scaling(cis_coverage,
                                                 only_fittable=True,
                                                 degree=degree,
                                                 n_cores=n_cores)
    rna_spline_scaling = _get_rna_scaling(cis_coverage,
                                          degree=degree,
                                          no_refine=no_refine,
                                          max_threshold=max_threshold,
                                          n_cores=n_cores)
    annotation = rdc_data.annotation
    refined_rna_splines = _refine_rna_splines(rna_spline_scaling, chrom_spline_scaling, annotation)
    rdc_data.write_scaling_batch(refined_rna_splines)
    rdc_data.is_scaling_fitted = True
    _rescale_rdc_data(rdc_data, fill_value, n_cores)
