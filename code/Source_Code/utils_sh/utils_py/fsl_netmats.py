import os
import sys
import argparse
import numpy as np
import pandas as pd
from fsl import nets
from datetime import datetime
from scipy.spatial.distance import squareform

def get_netmats(
        gen_outdir, ts_dir, parcel_type, 
        do_fulls=True, 
        do_partials=True, 
        do_amps=True, 
        debug=True,
        verbose=True
        ):
    start = datetime.now()

    subjfnames = [fname for fname in os.listdir(ts_dir) if "sub-" in fname]
    if subjfnames:
        sublist = [f.split('sub-')[1].split('.')[0] for f in subjfnames]
    else:
        # assumes only text files are in directory and all have form '<subjID>.txt'
        sublist = [f.split('.')[0] for f in os.listdir(ts_dir)]
    # load all subject data in given directory, sorted order
    ts = nets.load(ts_dir, 0.72, varnorm=0)

    if verbose:
        post_load = datetime.now()
        print(f"Elapsed time to load timeseries: {post_load - start}")
    
    if do_fulls:
        comp_netmats(
                ts, 
                gen_outdir,
                sublist,
                parcel_type,
                nm_type='F',
                verbose=verbose
                )

    if do_partials:
        comp_netmats(
                ts, 
                gen_outdir,
                sublist,
                parcel_type,
                nm_type='P',
                verbose=verbose
                )

    if do_amps:
        pre_comp = datetime.now()
        #try:
           #TS = np.squeeze(np.stack(ts.ts))
        #except ValueError:
            #return ts.ts
        subjfnames_debug = [fname for fname in os.listdir(ts_dir)]
        
        ts_data = []
        for fname in sorted(subjfnames_debug):
            ts_path = os.path.join(ts_dir, fname)
            ts_arr = np.loadtxt(ts_path)
            ts_data.append(ts_arr)
        
        amps_list = []
        for i, ts_arr in enumerate(ts_data):
            ts_arr = np.squeeze(ts_arr)
            amps = np.std(ts_arr, axis=0)
            amps_list.append(amps)
            if verbose:
                print(f"Amplitudes ({parcel_type}) dataset for subject {sublist[i]} has shape {amps.shape}")
        if debug:
            print(f"Number of subjects: {len(amps_list)}")
        if verbose:
            post_comp = datetime.now()
            print(f"Elapsed time to compute amplitudes: {post_comp - pre_comp}")
        
        #if debug:
            #print(f"Stack of loaded-in data has shape {TS.shape}")
        #if verbose:
            #post_comp = datetime.now()
            #print(f"Amplitudes ({parcel_type}) dataset has shape {amps.shape}")
            #print(f"Elapsed time to compute amplitudes: {post_comp - pre_comp}")

        #amps_list = [f"Amps_{parcel_type}_{i+1}" for i in range(amps.shape[1])]
        #amp_df = pd.DataFrame(data=amps, index=sublist, columns=amps_list)
        #amp_path = f"{gen_outdir}/Amplitudes_{parcel_type}.csv"
        amps_df = pd.DataFrame(data=amps_list, index=sublist)
        amp_path = os.path.join(gen_outdir, f"Amplitudes_{parcel_type}.csv")
        amps_df.to_csv(amp_path)


def comp_netmats(ts, gen_outdir, sublist, parcel_type, nm_type='F', verbose=True):
    pre_comp = datetime.now()
    if nm_type=='F':
        netmats = nets.netmats(ts, 'corr', True)
    elif nm_type=='P':
        netmats = nets.netmats(ts, 'ridgep', True, 0.1)
    else:
        raise ValueError("Invalid netmat computation type specified!")
    if verbose:
        post_comp = datetime.now()
        print(f"Elapsed time to compute partial netmats: {post_comp - pre_comp}")

    outdir = f"{gen_outdir}/{nm_type}netmats_{parcel_type}"
    save_subj_netmats(outdir, netmats, sublist)


def _make_square(vector, debug=False):
    n = np.sqrt(len(vector))
    _check_n(n)
    if n == int(n):

        if debug:
            ### debugging code ###
            print(f"reshaping input vector of length {n**2} into {int(n)}x{int(n)} matrix, error={n-int(n)}.")
            ### debugging code ###

        n = int(n)
        matrix = vector.reshape(n,n)
    else:
        matrix = squareform(vector, force='tomatrix', checks=True)
    return matrix

def _check_n(n):
    discriminant = np.sqrt(8*n**2 + 1)
    if discriminant % 2 == 1 and n > 1:
        n_triu = int((discriminant + 1)/2)
        print(f"AMBIGUOUS CASE: a vector of length {n**2} could be either a flattened {n}x{n} square matrix OR the upper-right triangle of a symmetric {n_triu}x{n_triu} matrix")

def save_subj_netmats(outdir, netmats, sublist, verbose=True):
    pre_save = datetime.now()

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    for i,sid in enumerate(sublist):
        ### !!! ###
        # assumes both sublist and netmats are sorted in subjID order!
        netvec = netmats[i]
        ### !!! ###
        netmat = _make_square(netvec)
        outpath = os.path.join(outdir, f"sub-{sid}.txt")
        np.savetxt(outpath, netmat)

    if verbose:
        post_save = datetime.now()
        print(f"Disk-written network matrices have shape: {_make_square(netmats[-1]).shape}")
        print(f"Elapsed time to write netmats to disk: {post_save - pre_save}")



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Binary classification: patient vs. control from resting-state data features")
    parser.add_argument(
            '-o',
            '--gen_outdir',
            type=str,
            default='/scratch/l.lexi/WAPIAW2024/Data/UKB/func_idp/partial_NMs',
            help='filepath to full set of IDPs (ie, indpendent variables) -- is directory in functional IDP case'
            )
    parser.add_argument(
            '-t',
            '--ts_dir',
            type=str,
            default='/scratch/l.lexi/WAPIAW2024/Data/UKB/func_idp/partial_NMs',
            help='filepath to full set of IDPs (ie, indpendent variables) -- is directory in functional IDP case'
            )
    parser.add_argument(
            '-p',
            '--parcel_type',
            type=str,
            default='/scratch/l.lexi/WAPIAW2024/Data/UKB/func_idp/partial_NMs',
            help='filepath to full set of IDPs (ie, indpendent variables) -- is directory in functional IDP case'
            )
    parser.add_argument(
            '-F',
            '--do_fulls',
            default=False,
            action='store_true',
            help="verbose flag: send problem parameters to output stream?"
            )
    parser.add_argument(
            '-P',
            '--do_partials',
            default=False,
            action='store_true',
            help="verbose flag: send problem parameters to output stream?"
            )
    parser.add_argument(
            '-A',
            '--do_amps',
            default=False,
            action='store_true',
            help="verbose flag: send problem parameters to output stream?"
            )
    args = parser.parse_args()

    get_netmats(
            args.gen_outdir, 
            args.ts_dir, 
            args.parcel_type, 
            do_fulls=args.do_fulls, 
            do_partials=args.do_partials, 
            do_amps=args.do_amps
            )

#test

### UPDATES:
# 1. faster output (i.e., not dataframe)
# 2. save out single-subject files (iterate)
# 3. save out each netmat as (symmetric) matrix
