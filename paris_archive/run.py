from agage_archive.config import Paths, data_file_path
from agage_archive.run import run_all
from pathlib import Path
from glob import glob
import os
import xarray as xr
import pandas as pd
import numpy as np

def preprocess():
    """Preprocess data files before running the main script    
    """
    current_dir = os.getcwd()

    # Sort ZSF first
    paths = Paths("paris", site="zsf")
    zsf_ms_folder = data_file_path("", "paris", sub_path=paths.gcms_flask_path)

    files = glob(str(zsf_ms_folder) + "/*_air.nc")
    for f in files:
        x = xr.open_dataset(f)
        # Subtract 92 min offset from Cedric Couret email to Joe Pitt 2025-06-07
        # 60 mins to put timestamp on UTC
        # 32 mins to get from the chromatogram time to the start of sampling
        x["sample_time"] = x["time"] - 92*60
        x_std = xr.open_dataset(f.split("_")[-2] + "_std.nc")
        spec_name = list(x.keys())[0].split("_")[0]
        rep = x_std[spec_name + "_C"].sel(time=slice(1721865600,4000000000)).std()  # average for latest std
        x[spec_name + "_std_stdev"] = x[spec_name + "_C"]
        x[spec_name + "_std_stdev"].values = np.repeat(float(rep), x.sizes["time"])
        x.to_netcdf(f + "_temp")
        os.system("rm -f " + f)
        os.system("mv " + f + "_temp " + f)

    # Now move ecd sf6 into the ms directory
    os.chdir(zsf_ms_folder)
    os.system("cp ../zugspitze-ecd/zsf_sf6.nc ./sf6_air.nc")
    x = xr.open_dataset("sf6_air.nc")
    x["SF6_C"] = x["sf6_C"]/1.002  # Convert to SIO-05 scale using Guillevic value
    x["SF6_std_stdev"] = x["sf6_C"]
    # Read in the repeatability values Cedric sent
    mf_rep = pd.read_csv("../zugspitze-ecd/SF6_Std_Stdv.txt",
                         sep="\t",
                         index_col="Date",
                         date_format="%d.%m.%Y")
    time_df = x["sf6_C"].to_dataframe()
    stdev_merge = pd.merge_ordered(time_df, mf_rep, fill_method="ffill", left_on="time", right_on="Date", how="left")["Stdv"]
    x["SF6_std_stdev"].values = stdev_merge
    x["time"] = x["time"].astype("int")// 10**9
    x["sample_time"] = x["time"] + 2500/2  # Need to add on 1/2 MS sample time
    x.to_netcdf("sf6_air.nc_temp")
    os.system("rm -f sf6_air.nc")
    os.system("mv sf6_air.nc_temp sf6_air.nc")

    os.chdir(current_dir)

    # Now sort BIR
    paths = Paths("paris", site="bir")
    bir_folder = data_file_path("", "paris", sub_path=paths.gcms_flask_path)

    files = glob(str(bir_folder) + "/*_air.nc")
    for f in files:
        x = xr.open_dataset(f)
        spec_name = list(x.keys())[5].split("_")[0]
        mean_stdev = x[spec_name + "_std_stdev"].mean(skipna=True)
        x[spec_name + "_std_stdev_orig"] = x[spec_name + "_std_stdev"].copy()
        # Earlier data has no std_stdev - use the mean of later data throughout
        x[spec_name + "_std_stdev"].values = np.repeat(float(mean_stdev), x.sizes["time"])
        x.to_netcdf(f + "_temp")
        os.system("rm -f " + f)
        os.system("mv " + f + "_temp " + f)

    # Now HUN
    paths = Paths("paris", site="hun")
    hun_folder = data_file_path("", "paris", sub_path=paths.gcms_flask_path)
    os.chdir(hun_folder)
    os.system("cp ../taunus-ecd_HUN_flask/sf6_air.nc ./sf6_air.nc")

    os.chdir(current_dir)

def postprocess():
    """Postprocess data files before running the main script 
    Need to convert filenames
    """

    # Change instrument_type (variable and attr) and comment for ZSF and CGR
    os.chdir("../data/paris")
    os.system("unzip paris-archive.zip -d paris-archive")
    for site in ["cgr", "zsf"]:
        files = glob("paris-archive/*/paris_" + site + "_*.nc")
        for f in files:
            x = xr.open_dataset(f)
            if site == "zsf" and x.attrs["species"] == "sf6":
                x["instrument_type"].values = np.repeat(6, x.sizes["time"])
                x["sampling_period"].values = np.repeat(3600, x.sizes["time"])
                x.attrs["instrument_type"] = "GCECD"
                x.attrs["comment"] = x.attrs["comment"].replace("GCMS Medusa flask", "GCECD")
                x.attrs["instrument"] = "Zugspitze GCECD"
                x.attrs["sampling_period"] = 3600
            else:
                x["instrument_type"].values = np.repeat(13, x.sizes["time"])
                x.attrs["instrument_type"] = "GCMS"
                x.attrs["comment"] = x.attrs["comment"].replace("GCMS Medusa flask", "GCMS")
            x.to_netcdf(f + "_temp")
            os.system("rm -f " + f)
            os.system("mv " + f + "_temp " + f)

    # Also change TOB SF6
    files = glob("paris-archive/sf6/paris_hun*.nc")
    x = xr.open_dataset(files[0])
    x["instrument_type"].values = np.repeat(-1, x.sizes["time"])
    x.attrs["instrument_type"] = "GCECD flask"
    x.attrs["comment"] = x.attrs["comment"].replace("GCMS Medusa flask", "GCECD flask")

    # Now rezip
    os.system("rm paris-archive.zip")
    os.chdir("paris-archive")
    os.system("zip -r ../paris-archive.zip .")

    # Remove ECD sf6 from MS dirs
    current_dir = os.getcwd()
    paths = Paths("paris", site="zsf")
    zsf_ms_folder = data_file_path("", "paris", sub_path=paths.gcms_flask_path)
    os.chdir(zsf_ms_folder)
    os.system("rm sf6_air.nc")

    os.chdir(current_dir)
    paths = Paths("paris", site="hun")
    hun_folder = data_file_path("", "paris", sub_path=paths.gcms_flask_path)
    os.chdir(hun_folder)
    os.system("rm sf6_air.nc")
    

if __name__ == "__main__":

    preprocess()

    run_all("paris",
            combined=False,
            baseline=False,
            monthly=False,
            top_level_only=True,
            resample = False)

    postprocess()
