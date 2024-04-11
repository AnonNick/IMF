import numpy as np
import calendar
import xarray as xr
from datetime import datetime, timedelta
import os
import ftplib


def check_if_made_of_ones(number):
            rounded_number = np.round(number, 2)  
            number_str = f"{rounded_number:.0f}"
            return all(char == '1' for char in number_str)


def is_invalid_sequence(day,datapoints):
    consecutive_invalid = 0
    count=0
    for point in datapoints:
        divided_by_9 = point / 9 
        all_nines = all(check_if_made_of_ones(num) for num in divided_by_9)
        if all_nines:
            consecutive_invalid += 1
            if consecutive_invalid >= 19:
                return True
        else:
            consecutive_invalid = 0
        count+=1
    return False

def get_data(year):
    backdate = 10
    data1 = np.loadtxt(f'/glade/work/nikhilr/tiegcm3.0/data/IMF/omni_asc/omni_min{year-1}.asc')
    data2 = np.loadtxt(f'/glade/work/nikhilr/tiegcm3.0/data/IMF/omni_asc/omni_min{year}.asc')
    data = np.concatenate((data1[-backdate:], data2), axis=0)
    

    last_valid_day = None

    unique_days = np.unique(data2[:,1])
    for day in unique_days:
        day_data = data2[data2[:,1] == day]  # Extract data for the day
        data_points = day_data[:,5:] 
        if is_invalid_sequence(day,data_points):
            last_valid_day = int(day) -1
            break
    if last_valid_day == None:
        last_valid_day = int(unique_days[-1])
    days = last_valid_day
    return data, days, backdate



def data_parse(data, year, days, backdate):    
    
    T = np.arange(60 * 24 * days) + 1  # Total time, minutes
    
    swden = data[:, 25]
    swvel = data[:, 21]
    swbx = data[:, 14]
    swby = data[:, 17]
    swbz = data[:, 18]

    timestamp_year = data[:,0]
    timestamp_day = data[:,1]
    timestamp_hour = data[:,2]
    timestamp_min = data[:,3]

    t = data[:, 26]

    swden[swden > 500] = np.nan
    swvel[swvel > 1e4] = np.nan
    swbx[swbx > 1e3] = np.nan
    swby[swby > 1e3] = np.nan
    swbz[swbz > 1e3] = np.nan

    timestamp = [''] * len(T)
    swbx4 = np.zeros(len(T))
    swbx4[:] = np.nan
    swby4 = np.zeros(len(T))
    swby4[:] = np.nan
    swbz4 = np.zeros(len(T))
    swbz4[:] = np.nan
    swden4 = np.zeros(len(T))
    swden4[:] = np.nan
    swvel4 = np.zeros(len(T))
    swvel4[:] = np.nan
    velMask = np.ones(len(T))  
    denMask = np.ones(len(T))  
    bxMask = np.ones(len(T))   
    byMask = np.ones(len(T))   
    bzMask = np.ones(len(T))   
    date = np.zeros(len(T))

    count = 0  
    print("    Timestamps")
    for i in T:

        ts_year = int(timestamp_year[count + backdate])
        day_of_year = int(timestamp_day[count + backdate])
        hour = int(timestamp_hour[count + backdate])
        minute = int(timestamp_min[count + backdate])
        timestamp[count]=(datetime(ts_year, 1, 1) + timedelta(days=day_of_year - 1, hours=hour, minutes=minute)).strftime('%Y-%m-%dT%H:%M:%S')
        if count%(60) == 0:
            print(f'        {timestamp[count]}')
        Dstart = i -1
        Dend = i + backdate - 1
        if not np.isnan(swbx[Dstart:Dend]).all():  # Checks if not all values are NaN
            swbx4[count] = np.nanmean(swbx[Dstart:Dend])
        else:
            swbx4[count] = np.nan  # Or any other default value as required
        # Repeat the check for other variables (swby, swbz, swden, swvel)
        if not np.isnan(swby[Dstart:Dend]).all():
            swby4[count] = np.nanmean(swby[Dstart:Dend])
        else:
            swby4[count] = np.nan  # Default value for all-NaN slices

        if not np.isnan(swbz[Dstart:Dend]).all():
            swbz4[count] = np.nanmean(swbz[Dstart:Dend])
        else:
            swbz4[count] = np.nan  # Default value for all-NaN slices

        if not np.isnan(swden[Dstart:Dend]).all():
            swden4[count] = np.nanmean(swden[Dstart:Dend])
        else:
            swden4[count] = np.nan  # Default value for all-NaN slices

        if not np.isnan(swvel[Dstart:Dend]).all():
            swvel4[count] = np.nanmean(swvel[Dstart:Dend])
        else:
            swvel4[count] = np.nan  # Default value for all-NaN slices
        date[count] = year*1000 + 1 + (count)/(60*24)
        count += 1


    swbx4_nan_indices = np.where(np.isnan(swbx4))[0]
    bxMask[swbx4_nan_indices] = 0
    swbx4nan = swbx4[~np.isnan(swbx4)]
    Tbx4 = T[~np.isnan(swbx4)]
    swbx41 = np.interp(T, Tbx4, swbx4nan)
    swbx4 = swbx41
    
    swby4_nan_indices = np.where(np.isnan(swby4))[0]
    byMask[swby4_nan_indices] = 0
    swby4nan = swby4[~np.isnan(swby4)]
    Tby4 = T[~np.isnan(swby4)]
    swby41 = np.interp(T, Tby4, swby4nan)
    swby4 = swby41
    
    swbz4_nan_indices = np.where(np.isnan(swbz4))[0]
    bzMask[swbz4_nan_indices] = 0
    swbz4nan = swbz4[~np.isnan(swbz4)]
    Tbz4 = T[~np.isnan(swbz4)]
    swbz41 = np.interp(T, Tbz4, swbz4nan)
    swbz4 = swbz41
    
    swden4_nan_indices = np.where(np.isnan(swden4))[0]
    denMask[swden4_nan_indices] = 0
    swden4nan = swden4[~np.isnan(swden4)]
    Tden4 = T[~np.isnan(swden4)]
    swden41 = np.interp(T, Tden4, swden4nan)
    swden4 = swden41
    
    swvel4_nan_indices = np.where(np.isnan(swvel4))[0]
    velMask[swvel4_nan_indices] = 0
    swvel4nan = swvel4[~np.isnan(swvel4)]
    Tvel4 = T[~np.isnan(swvel4)]
    swvel41 = np.interp(T, Tvel4, swvel4nan)
    swvel4 = swvel41

    return swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp
   
def save_netcdf (file_path,swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp):
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
    ds.attrs["Description"] = "10-minute average of OMNI data trailed by 1 minutes. Sampled to minute output"
    ds.attrs["CreationTime"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f") 
    ds.attrs["Source"] = "Hourly OMNI combined 1AU IP Data"
    ds.attrs["Version"] = "1.0.0"
    ds.attrs["CreatedBy"] = "nikhilr"  
    ds.attrs["url_reference"] = "https://omniweb.gsfc.nasa.gov/ow_min.html"


    ds.to_netcdf(file_path)

    print(f"    Dataset has been saved to {file_path}")

def should_download(file_mod_time, current_time):
    # Check if we are past the 13th of the current month
    if current_time.day > 13:
        threshold_time = current_time.replace(day=13)
    else:
        last_month = current_time.replace(day=1) - timedelta(days=1)
        threshold_time = last_month.replace(day=20)
    
    return file_mod_time > threshold_time    

def download_files():
    ftp_url = 'spdf.gsfc.nasa.gov'
    directory = '/pub/data/omni/high_res_omni/'
    local_directory = '/glade/work/nikhilr/tiegcm3.0/data/IMF/omni_asc'
    os.makedirs(local_directory, exist_ok=True)
    
    with ftplib.FTP_TLS(ftp_url) as ftp:
        ftp.login() 
        ftp.prot_p() 
        ftp.cwd(directory) 
        
        files_to_download = ftp.nlst('omni_min*.asc')
        
        current_time = datetime.utcnow()
        years_downloaded = []

        for filename in files_to_download:
            local_path = os.path.join(local_directory, filename)
            file_exists_locally = os.path.exists(local_path)
            year = filename[-8:-4]
            if not file_exists_locally:
                # If file does not exist locally, download it without checking modification time
                with open(local_path, 'wb') as local_file:
                    ftp.retrbinary(f'RETR {filename}', local_file.write)
                years_downloaded.append(int(year))
            else:
                # If file exists, check modification time
                mod_time_str = ftp.sendcmd(f'MDTM {filename}')[4:].strip()
                file_mod_time = datetime.strptime(mod_time_str, '%Y%m%d%H%M%S')

                if should_download(file_mod_time, current_time):
                    with open(local_path, 'wb') as local_file:
                        ftp.retrbinary(f'RETR {filename}', local_file.write)
                    years_downloaded.append(int(year))
    
    return years_downloaded



def main():
    years= download_files()
    year_limit = 1982
    years = [year for year in years if year >= year_limit]
    print(f"Updating {years}")
    #years = list(range(1982, 2000 + 1))
    #years = [2024]
    for year in years:
        print(f"Year: {year}")
        data, days, backdate=get_data(year)
        swbx4, bxMask, swby4, byMask, swbz4, bzMask, swden4, denMask, swvel4, velMask, date, timestamp = data_parse(data, year, days, backdate)
        file_path = '/glade/work/nikhilr/tiegcm3.0/data/IMF'
        file_name = f'imf_OMNI_{year}001-{year}{str(days).zfill(3)}.nc'
        file_path = file_path+'/'+file_name
        save_netcdf (file_path,swbx4,bxMask,swby4,byMask,swbz4,bzMask,swden4,denMask,swvel4,velMask,date,timestamp)
    
if __name__ == "__main__":
    main()