# This function is used to modify the date for a different year than the NEI emissions file
# Created by Madankui Tao on July 24, 2023

# Input: Hourly mean NEI filename (in the format of '')
# Output: NEI filename for the target year

#=====Dependent library & functions======
#--------------------------------------------------------------------------
import os
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from datetime import datetime
def extract_date_and_time(filename):
    # Extract the date and time part of the filename
    date_time_str = filename[13:]
    # Convert the date and time string to a datetime object (without seconds)
    date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d_%H:%M:%S')
    return date_time_obj
# # Sort the list based on date and time using the custom key function
# sorted_filenames = sorted(filenames, key=extract_date_and_time)

#=====The function======
#--------------------------------------------------------------------------
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

def NEIfileDatetime_shift_BYweekofday(inpNEIfile, outyear):
    """
    Adjusts the date in a NEI file to match the same day in a week in the specified output year.

    Parameters:
        inpNEIfile (str): The NEI file with the input datetime in the format 'wrfchemi_d01_YYYY-MM-DD_HH:MM:SS'.
        outyear (int): The target output year for the adjusted datetime.

    Returns:
        str: The adjusted NEI filename for the specified output year in the format 'wrfchemi_d01_YYYY-MM-DD_HH:MM:SS'.
    """
    # Takes NEI file for the input datetime
    inp_date_time_obj = datetime.strptime(inpNEIfile[13:], '%Y-%m-%d_%H:%M:%S')

    # get the years between input and output
    nyears = outyear - inp_date_time_obj.year

    # find the date in the target year
    outp_date_time_obj = inp_date_time_obj + relativedelta(years=nyears)

    # adjusted for week in a day
    delta_days = outp_date_time_obj - inp_date_time_obj

    # Calculate the number of weeks to shift the weekday
    num_weeks_shift = delta_days.days // 7

    # Calculate the number of remaining days after removing whole weeks
    remaining_days = delta_days.days % 7

    # then will shift the date by 'remaining_days' to match the same day in a week
    adjt_outp_date_time_obj = outp_date_time_obj + timedelta(days=-remaining_days)

    # Format the datetime object to get the full weekday names of the original and adjusted dates
    original_weekday = inp_date_time_obj.strftime('%A')
    adjusted_weekday = adjt_outp_date_time_obj.strftime('%A')

    # Print the day of the week of the original and adjusted dates for verification
    print(f"{original_weekday} matched with {adjusted_weekday}")

    # return the NEI filename after adjustment
    outpNEIfile = 'wrfchemi_d01_' + adjt_outp_date_time_obj.strftime('%Y-%m-%d_%H:%M:%S')

    return outpNEIfile

# #=====Example to use this function======
# #--------------------------------------------------------------------------
# inpNEIfile = 'wrfchemi_d01_2017-08-01_01:00:00'
# outyear = 2020

# outpNEIfile = NEIfileDatetime_shift_BYweekofday(inpNEIfile, outyear)

# #=====Function to locate the 2017 NEI filename given the desired output date======
# #--------------------------------------------------------------------------

def NEIfileDatetime_shift_BACKWARDto2017(outpNEIfile_datetime, outyear):
    """
    Calculate the original 2017 datetime for a given NEI file with adjusted datetime and the specified output year.

    Parameters:
        outpNEIfile_datetime (str): The NEI file with the adjusted datetime in the format 'YYYY-MM-DD_HH:MM:SS'.
        outyear (int): The target output year for the adjusted datetime.

    Returns:
        str: The original NEI filename datetime for the specified output year in the format 'YYYY-MM-DD_HH:MM:SS'.
    """
    # Takes outpNEIfile_datetime for the adjusted datetime
    outp_date_time_obj = datetime.strptime(outpNEIfile_datetime, '%Y-%m-%d_%H:%M:%S')

    # Calculate the number of years between the output year and 2017
    nyears = outp_date_time_obj.year - 2017

    # Find the date in the input year (2017)
    inp_date_time_obj = outp_date_time_obj - relativedelta(years=nyears)

    # Adjust for the week in a day
    delta_days = outp_date_time_obj - inp_date_time_obj

    # Calculate the number of weeks to shift the weekday
    num_weeks_shift = delta_days.days // 7

    # Calculate the number of remaining days after removing whole weeks
    remaining_days = delta_days.days % 7

    # Shift the date by 'remaining_days' to match the same day in a week
    adjt_inp_date_time_obj = inp_date_time_obj + timedelta(days=remaining_days)

    # Format the datetime object to get the full weekday names of the original and adjusted dates
    original_weekday = adjt_inp_date_time_obj.strftime('%A')
    adjusted_weekday = outp_date_time_obj.strftime('%A')

    # Print the day of the week of the original and adjusted dates for verification
    print(f"{adjusted_weekday} matched with {original_weekday}")

    # Return the NEI filename after backward adjustment
    inpNEIfile_datetime = adjt_inp_date_time_obj.strftime('%Y-%m-%d_%H:%M:%S')

    return inpNEIfile_datetime
