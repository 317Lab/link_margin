import numpy as np
import pandas as pd
def read_gps_file(filename):
    df = pd.read_excel(filename)
    # strip white space from column names
    df.columns = df.columns.str.strip()
    df = df[['Flight Time', 'Lat', 'Long', 'Alt']]
    df[['Alt']]=df[['Alt']]*1000
    arr = df.to_numpy()
    # last three rows have negative time, removed for safety
    arr = np.delete(arr, np.s_[len(arr)-3:len(arr)], axis=0)
    # shift times to be positive
    arr[:,0] = arr[:,0] + np.abs(arr[0,0])

    times = arr[:,0]
    lla = arr[:,1:arr.shape[1]]
    return times, lla

file = "giraff_381_sub_gps.xlsx"
savepath = "381_sub_processed"
times, lla = read_gps_file(file)
np.savez_compressed(savepath, times=times, lla=lla)