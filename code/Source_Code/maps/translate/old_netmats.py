import os
import argparse
import numpy as np
import pandas as pd
from itertools import chain
from scipy.linalg import block_diag

# produces a block matrix from map
def make_fIDP_Yeo_mask(thresh, map_dir, parcel_type="Schaefer"):
    maskpath = os.path.join(map_dir, f"{parcel_type}_Yeo_mask.csv")

    # assumes [functional parcellation]-->Yeo map is saved according to strict naming convention
    map_df = pd.read_csv(os.path.join(map_dir, f"{parcel_type}_map_Yeo.csv"), index_col=0)
    mask = _get_mask(thresh, map_df)
    return mask 

def _get_mask(thresh, map_df, debug=False, verbose=False):
    thresh_df = map_df[map_df["Overlap_percent"] > thresh]

    blkdiag_list = []
    ordered_index = []
    for parcel_num in thresh_df["Yeo_parcels"].unique():
        parcel_df = thresh_df[thresh_df["Yeo_parcels"] == parcel_num]

        parcel_idx = list(parcel_df.index)
        ordered_index.append(parcel_idx)
        
        n_parcel = len(parcel_idx)
        blkdiag_list.append(parcel_num*np.ones((n_parcel,n_parcel)))

    unordered_block = block_diag(*blkdiag_list)
    if verbose:
        print(f"index subgroups have lengths: {[len(i) for i in ordered_index]}")
        print(f"matrix blocks have shapes: {[i.shape for i in blkdiag_list]}")
        print(f"the (unordered) block diagonal matrix has shape: {unordered_block.shape}")

    idx = [i-1 for i in list(chain.from_iterable(ordered_index))]   # parcel indices start at 1, python indices start from 0 
    mask = _replace_submatrix(np.zeros((len(map_df),len(map_df))), idx, idx, unordered_block)
    
    if debug:
        ### debugging code ###
        print(f"mask has dimensions: {mask.shape}")
        return blkdiag_list, unordered_block, map_df, ordered_index
        ### debugging code ###

    np.fill_diagonal(mask, 0)
    
    return mask

# replace a submatrix (specified by ind1, ind2) of given matrix "mat" with "mat_replace"
def _replace_submatrix(mat, ind1, ind2, mat_replace):
    for i, index in enumerate(ind1):
        mat[index, ind2] = mat_replace[i, :]
        
    return mat

# for given parcellation map, estimates output as function of overlap threshold
def _thresh_loss(df, thresh, parcel_type="Schaefer"):
    parcelprop = 1 - len(df[df["Overlap_percent"] < thresh])/len(df)
    netmatprop = parcelprop**2
    print(f"In the {parcel_type} parcellation, requiring at least {thresh}% Yeo overlap retains:")
    print(f"{np.round(100*parcelprop,1)}% of parcels")
    print(f"{np.round(100*netmatprop,1)}% of netmat components")


# save mask out given output directory and parcellation name type
def save_mask(mask, outdir, thresh=50, parcel_type="Schaefer", debug=True):
    mask_outname = f"mask{thresh}_{parcel_type}_to_Yeo.txt"
    flat_outname = f"halfmask{thresh}nonzero_{parcel_type}_to_Yeo.txt"

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    if debug:
        ### debugging code ###
        print(f"mask has dimensions: {mask.shape}")
        ### debugging code ###

    np.savetxt(os.path.join(outdir, mask_outname), mask, fmt='%i')

    flat_mask = _flatmask_from_sym(mask, mask)
    np.savetxt(os.path.join(outdir, flat_outname), flat_mask, fmt='%i')

# given a symmetric matrix A, pull the entries of its upper right triangle corresponding to non-zero entries in mask
def _flatmask_from_sym(A, mask, keep_diags=False):
    assert A.shape == mask.shape, "Input symmetric matrix and mask must have same shape"
    if keep_diags:
        k=0
    else:
        k=1

    if len(A.shape) > 1:
        assert np.allclose(A,A.T), "If input is a matrix, it must be symmetric (and, implicitly, square)"
        A = A[np.triu_indices(A.shape[0],k)]
        mask = mask[np.triu_indices(mask.shape[0],k)]

    A_condensed = A[np.nonzero(mask)]
    return A_condensed


def mask_dataset(indir, thresh, debug=True):
    if "Glasser" in indir:
        parcel_type="Glasser"
    elif "Schaefer" in indir:
        parcel_type="Schaefer"
    else:
        raise ValueError("Input directory must contain substring specifying parcel type")

    # name the output directory (whatever you wanna do, i just put something here)
    outdir = indir + f"flat{thresh}_yeo"

    # check for mask or make it or whatever
    mask = ??? something using make_fIDP_yeo_mask(???, thresh, parcel_type=parcel_type)

    # apply mask on per-subject basis:
    for datafile in os.listdir(indir):
        netmat = ??? some_load_function(os.path.join(indir,datafile))   # not sure if its best as pd.read_csv, np.genfromtxt (and/or with what delimiter), should check data
        netmat_condensed = _flatmask_from_sym(netmat, mask)
        some_save_function(os.path.join(outdir, datafile))              # same thoughts as the comment for loading stuff in


## PARSER
##################################################################################################################################
# parses inputs, generates derivative intermediates, loads and splits data, fits and predicts with RF classifiers, and saves results
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Binary classification: patient vs. control from resting-state data features")
    parser.add_argument(
            '-t',
            '--overlap_thresh',
            type=float,
            default=50,
            help='minimal allowable percentage of overlap between given parcel and its 1st-place yeo component match'
            )
    parser.add_argument(
            '-d',
            '--map_dir',
            type=str,
            default="/scratch/l.lexi/WAPIAW2024/Source_Code/maps",
            help='directory containing maps to Yeo parcellation (from other parcellations)'
            )
    parser.add_argument(
            '-o',
            '--output_dir',
            type=str,
            default="/scratch/l.lexi/WAPIAW2024/Source_Code/maps/netmat_masks",
            help='directory in which to output'
            )

    args = parser.parse_args()

    ptypes = ["Glasser", "Schaefer"]
    for parcel in ptypes:
        mask = make_fIDP_Yeo_mask(
                args.overlap_thresh, 
                args.map_dir,
                parcel_type = parcel
                )

        save_mask(
                mask,
                args.output_dir,
                thresh = args.overlap_thresh,
                parcel_type = parcel
                )
    #code to get half flat mask for folder
    #change these for each run
    #    indir = "/scratch/l.lexi/WAPIAW2024/Data/UKB/func_idp/Fnetmats_glasser"

####   See 'def mask_dataset(indir, thresh):' (line 101) for pseudocode that generalizes the steps below   ####

#   indir = "/scratch/l.lexi/WAPIAW2024/Data/UKB/func_idp/archived/NetMats"
#   thresh = args.overlap_thresh
#   outdir = indir + f"_flat_yeo{thresh}"
#   for nm in os.listdir(indir):
#       map_df = pd.read_csv(os.path.join(indir, nm), index_col=0)
#       mask = _get_mask(args.overlap_thresh, map_df)
#       half_mask = np.triu(mask)
#       flat_half_mask = half_mask.flatten()
#       flattened_non_zero = flat_half_mask[flat_half_mask != 0]
#       # save flat mask
#       nm_name = nm.split(".")[0]
#       mask_outname = f"{nm_name}_flat_yeo{thresh}.txt"
#       if not os.path.isdir(outdir):
#               os.mkdir(outdir)
#       outpath = os.path.join(outdir, mask_outname)
#       np.savetxt(outpath, flattened_non_zero, fmt='%i')
