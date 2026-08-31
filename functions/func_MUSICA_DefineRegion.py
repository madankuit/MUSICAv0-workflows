"""
func_MUSICA_DefineRegion.py
Defines geographic region boundaries and timezone offsets for model evaluation
over the contiguous United States (CONUS) and sub-regions.

Provides:
    - latlonbound(fileregion): lat/lon boundaries for named US regions
    - timezone_approx1330LT(fileregion): UTC time window for approx. 1:30 PM local time
    - UTCtimezone_offset_summer(fileregion): UTC offset (positive hours) for summer/DST
    - regional_UTCtimezone_offset_summer(fileregion): UTC offset (negative hours) for summer/DST

Supported regions: WestCoast, Mountain, Midwest, Southwest, Southeast, Northeast, CONUS

MODIFICATION HISTORY:
    Initial version
"""

# Define the region boundaries

#=====lat.lon bounds======
#--------------------------------------------------------------------------

# CONUS LAT LONG bound (based on NEI boundary)
# # 40 deg LAT by 80 deg LONG
# CONUS_LAT_UP = 60
# CONUS_LAT_BOT = 15
# CONUS_LONG_LEFT = -140
# CONUS_LONG_RIGHT = -50
### Past

def latlonbound(fileregion):
    """This function returns the lats and lons boundary of a given region
        Regions include: 'WestCoast', 'Mountain', 'Midwest', 'Southwest', 'Southeast', 'Northeast'
    """
    # Define regions        
    if fileregion == "WestCoast":
        lon_right,lat_bot,lon_left,lat_up = [-125, 31, -115, 49]
    elif fileregion == "Mountain":
        lon_right,lat_bot,lon_left,lat_up = [-115, 31, -103, 49]
    elif fileregion == "Midwest":
        lon_right,lat_bot,lon_left,lat_up = [-103, 36.5, -85, 49]
    elif fileregion == "Southwest":
        lon_right,lat_bot,lon_left,lat_up = [-103, 25.5, -85, 36.5]
    elif fileregion == "Southeast":
        lon_right,lat_bot,lon_left,lat_up = [-85, 24.5, -75, 36.5]
    elif fileregion == "Northeast":
        lon_right,lat_bot,lon_left,lat_up = [-85, 36.5, -67, 49]
    elif fileregion == "CONUS":
        lon_right,lat_bot,lon_left,lat_up = [-126, 23, -66, 50]
    
    return [lon_right,lat_bot,lon_left,lat_up]

## Divide the U.S.
CONUS_LONG_LEFT, CONUS_LAT_BOT, CONUS_LONG_RIGHT, CONUS_LAT_UP = latlonbound('CONUS')

#=====Region labels======
#--------------------------------------------------------------------------

# West Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("WestCoast")
WestUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
WestUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelWestCoast = 'West Coast ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'
# print(labelWestCoast) # uncomment if need to verify if the info has been loaded

# Mountain Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("Mountain")
MountainUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
MountainUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelMountain = 'Mountain ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'

# Midwest Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("Midwest")
MidwestUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
MidwestUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelMidwest = 'Midwest ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'

# Southwest Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("Southwest")
SouthwestUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
SouthwestUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelSouthwest = 'Southwest ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'

# Southeast Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("Southeast")
SoutheastUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
SoutheastUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelSoutheast = 'Southeast ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'

# Northeast Region
lon_right,lat_bot,lon_left,lat_up = latlonbound("Northeast")
NortheastUS_lon = [lon_right, lon_right, lon_left, lon_left, lon_right]
NortheastUS_lat = [lat_bot, lat_up, lat_up, lat_bot, lat_bot]
labelNortheast = 'Northeast ('+str(lat_bot)+'-'+str(lat_up)+'°N, '+str(-lon_left)+'-'+str(-lon_right)+'°W)'

#=====Start and end time======
#--------------------------------------------------------------------------
# Specify time for the given region
def timezone_approx1330LT(fileregion):
    """This function returns the approximate time selection based on the timezone of a given region
       *With daylight saving  
       Select for 1 and 2 PM Local Time
       Options: WestCoast,Mountain,Midwest,Southwest,Southeast,Northeast
    """
    # Define regions        
    if fileregion == "WestCoast":
        # Pacific (UTC-7)
        averageTimeStart,averageTimeEnd = ['T20:00:00', 'T21:00:00']
    elif fileregion == "Mountain":
        # Mountain (UTC-6)
        averageTimeStart,averageTimeEnd = ['T19:00:00', 'T20:00:00']
    elif fileregion == "Midwest":
        # Central (UTC-5)
        averageTimeStart,averageTimeEnd = ['T18:00:00', 'T19:00:00']
    elif fileregion == "Southwest":
        # Central (UTC-5)
        averageTimeStart,averageTimeEnd = ['T18:00:00', 'T19:00:00']
    elif fileregion == "Southeast":
        # Eastern (UTC-4)
        averageTimeStart,averageTimeEnd = ['T17:00:00', 'T18:00:00']
    elif fileregion == "Northeast":
        # Eastern (UTC-4)
        averageTimeStart,averageTimeEnd = ['T17:00:00', 'T18:00:00']
    
    return [averageTimeStart,averageTimeEnd]

#--------------------------------------------------------------------------
# Specify timezone for the given region
def UTCtimezone_offset_summer(fileregion):
    """This function returns the timezone of a given region
       *With daylight saving  
       Options: WestCoast,Mountain,Midwest,Southwest,Southeast,Northeast
    """
    # Define regions        
    if fileregion == "WestCoast":
        # Pacific (UTC-7)
        UTCtimezone_offset = 7
    elif fileregion == "Mountain":
        # Mountain (UTC-6)
        UTCtimezone_offset = 6
    elif fileregion == "Midwest":
        # Central (UTC-5)
        UTCtimezone_offset = 5
    elif fileregion == "Southwest":
        # Central (UTC-5)
        UTCtimezone_offset = 5
    elif fileregion == "Southeast":
        # Eastern (UTC-4)
        UTCtimezone_offset = 4
    elif fileregion == "Northeast":
        # Eastern (UTC-4)
        UTCtimezone_offset = 4
    
    return UTCtimezone_offset

# Specify timezone for the given region
def regional_UTCtimezone_offset_summer(fileregion):
    """This function returns the timezone of a given region
       *With daylight saving  
       Options: WestCoast,Mountain,Midwest,Southwest,Southeast,Northeast
    """
    # Define regions        
    if fileregion == "WestCoast":
        # Pacific (UTC-7)
        UTCtimezone_offset = -7
    elif fileregion == "Mountain":
        # Mountain (UTC-6)
        UTCtimezone_offset = -6
    elif fileregion == "Midwest":
        # Central (UTC-5)
        UTCtimezone_offset = -5
    elif fileregion == "Southwest":
        # Central (UTC-5)
        UTCtimezone_offset = -5
    elif fileregion == "Southeast":
        # Eastern (UTC-4)
        UTCtimezone_offset = -4
    elif fileregion == "Northeast":
        # Eastern (UTC-4)
        UTCtimezone_offset = -4
    
    return UTCtimezone_offset