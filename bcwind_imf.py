import h5py
import numpy as np
import xarray as xr
import pandas as pd
import argparse
from datetime import datetime, timedelta
def read_hdf5(file_path):
    data = {}
    with h5py.File(file_path, 'r') as f:
        for key in f.keys():
            data[key] = f[key][()]
    return data

def process_data(data):
    swbx4 = data['Bx']
    swby4 = data['By']
    swbz4 = data['Bz']
    swden4 = data['D']
    swvel4 = data['Va']
    ut_str = [x.decode('utf-8') for x in data['UT']]
    timestamp = pd.to_datetime(ut_str)

    # Example processing, adjust as needed
    swden4[swden4 > 500] = np.nan
    swvel4[swvel4 > 1e4] = np.nan
    swbx4[swbx4 > 1e4] = np.nan
    swby4[swby4 > 1e4] = np.nan
    swbz4[swbz4 > 1e4] = np.nan
    
    velMask = np.ones(len(timestamp))  
    denMask = np.ones(len(timestamp))  
    bxMask = np.ones(len(timestamp))   
    byMask = np.ones(len(timestamp))   
    bzMask = np.ones(len(timestamp))   
    date = np.zeros(len(timestamp))

    count = 0
    firstday = timestamp[0].dayofyear
    daycount = ((firstday-1 )*(60 * 24)) + timestamp[0].hour*60 + timestamp[0].minute
    for time in timestamp:
        year = time.year
        day_of_year = time.dayofyear
        hour = time.hour
        minute = time.minute
        #date[count] = round(year*1000 + day_of_year + (hour*60+minute)/(60*24),8)
        date[count] = round(year*1000 + 1 + (daycount)/(60*24),8)
        print (time,date[count])
        count = count + 1
        daycount = daycount + 1
    
    return swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp


def save_netcdf (hdf5_file_path,file_path,swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp):
    ndata = len(date)  
    ds = xr.Dataset(
        {
            "bx": ("ndata", swbx4, {"units": "nT", "long_name": "IMF Bx"}),
            "bxMask": ("ndata", bxMask.astype('int8'), {"units": "boolean", "long_name": "Quality flag: 0=data derived from linear interpolation."}),
            "by": ("ndata", swby4, {"units": "nT", "long_name": "IMF By"}),
            "byMask": ("ndata", byMask.astype('int8'), {"units": "boolean", "long_name": "Quality flag: 0=data derived from linear interpolation."}),
            "bz": ("ndata", swbz4, {"units": "nT", "long_name": "IMF Bz"}),
            "bzMask": ("ndata", bzMask.astype('int8'), {"units": "boolean", "long_name": "Quality flag: 0=data derived from linear interpolation."}),
            "swden": ("ndata", swden4, {"units": "cm^{-3}", "long_name": "solar wind density"}),
            "denMask": ("ndata", denMask.astype('int8'), {"units": "boolean", "long_name": "Quality flag: 0=data derived from linear interpolation."}),
            "swvel": ("ndata", swvel4, {"units": "km/s", "long_name": "solar wind velocity"}),
            "velMask": ("ndata", velMask.astype('int8'), {"units": "boolean", "long_name": "Quality flag: 0=data derived from linear interpolation."}),
            "date": ("ndata", date, {"long_name": "year-day plus fractional day: yyyyddd.frac"}),
            "timestamp": ("ndata", timestamp, {"long_name": "Timestamp of the data: YYYY-MM-DDTHH:MM:SS"}),
        },
        coords={"ndata": np.arange(ndata)}
    )

    # Global attributes
    ds.attrs["Description"] = "BCWIND.h5 to minute output IMF data"
    ds.attrs["CreationTime"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f") 
    ds.attrs["Source"] = f"{hdf5_file_path}"
    ds.attrs["Version"] = "1.0.0"
    ds.attrs["CreatedBy"] = "nikhilr"  
    ds.attrs["url_reference"] = "https://github.com/AnonNick/IMF"


    ds.to_netcdf(file_path)

    print(f"    Dataset has been saved to {file_path}")



def main():
    parser = argparse.ArgumentParser(description="Convert bcwind HDF5 file to IMF NetCDF file.")
    parser.add_argument('hdf5_file_path', type=str, help="Path to the bcwind HDF5 file")
    args = parser.parse_args()
    hdf5_file_path = args.hdf5_file_path
    data = read_hdf5(hdf5_file_path)
    swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp = process_data(data)
    
    from_date = timestamp.min().strftime('%Y%j')
    to_date = timestamp.max().strftime('%Y%j')
    
    output_file_path = f'imf_bcwind_{from_date}-{to_date}.nc'  # Construct the output filename
    save_netcdf(hdf5_file_path,output_file_path, swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp)
    print(f"File saved as {output_file_path}")

if __name__ == "__main__":
    main()
