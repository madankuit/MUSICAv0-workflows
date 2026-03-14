'''
ModelEval_statistical_tests.py
this script contains functions of statistical tests designed for model evaluation, i.e., comparing model simulations with SLAMS or TROPOMI
MODIFICATION HISTORY:
    Madankui Tao, 3, OCT, 2023: VERSION 1.0
    - Initial version
    Madankui Tao, 17, DEC, 2023: VERSION 1.1
    - Add options to specify p-value and print number format
'''
### Module import ###
import numpy as np # for array manipulation and basic scientific calculation
import xarray as xr # To read NetCDF files

#================================================================================================
### R2 for RMA Regression ###
from pylr2 import regress2
def ReducedMajorAxisRegression_R2_TwoArComp(data1, data2):
    """Returns r2 of two arrays for comparison"""
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]

    results = regress2(nonNA_data1, nonNA_data2, _method_type_2="reduced major axis")
    R2 = np.square(results['r'])
    
    print(r"R2 = %5.2f" % (R2))
    
    # Meandiff = np.nanmean(nonNA_data1-nonNA_data2)
    # Meandiff_percent = np.nanmean(((nonNA_data1-nonNA_data2)/nonNA_data1)*100)
    # print(r"Mean Diff Data1 - 2 = %5.2f; %5.2f percent" % (Meandiff,Meandiff_percent))
    
    return R2

#================================================================================================
### greater/less two datasets ###
from scipy import stats
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html#scipy.stats.wilcoxon
def Oneside_Wilcoxon_test_result(Datestr1, data1, Datestr2, data2, givenpval, numbertype):
    """ T-test on two datasets for: 
        less and greater
        givenpval: the p-value threshold
        numbertype: the format of printed results
    """
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]
    
    df = len(nonNA_data1)+len(nonNA_data2)-2
    print("Degree of freedom =",df)
    # the distribution underlying d is stochastically less than a distribution symmetric about zero.
    sttsval, pvalsless = stats.wilcoxon(nonNA_data1, nonNA_data2, alternative='less') #"two-sided"
    print("%s < %s, p-value = %5.4f" % (Datestr1,Datestr2,np.nanmean(pvalsless)))
    # the distribution underlying d is stochastically greater than a distribution symmetric about zero
    sttsval, pvalsgreater = stats.wilcoxon(nonNA_data1, nonNA_data2, alternative='greater') #"two-sided"
    print("%s > %s, p-value = %5.4f" % (Datestr1,Datestr2,np.nanmean(pvalsgreater)))
    if pvalsless<givenpval or pvalsgreater<givenpval:
        print("The pseudomedian is NOT located at zero.")
    
    # calculate the difference
    Meandiff = np.nanmean(nonNA_data1-nonNA_data2)
    Meandiff_percent = Meandiff/np.nanmean(nonNA_data1)*100
    Mediandiff = np.nanmedian(nonNA_data1-nonNA_data2)
    
    # with specified format, print the results
    if numbertype=='f':
        print(r"Mean Diff %s - %s = %5.2f; %5.0f percent" % (Datestr1,Datestr2,Meandiff,Meandiff_percent))
        print(r"Median of Diff %s - %s = %5.2f" % (Datestr1,Datestr2,Mediandiff))
    elif numbertype=='e':
        print(r"Mean Diff %s - %s = %5.2e; %5.0f percent" % (Datestr1,Datestr2,Meandiff,Meandiff_percent))
        print(r"Median of Diff %s - %s = %5.2e" % (Datestr1,Datestr2,Mediandiff))
    
    return Meandiff,Meandiff_percent,Mediandiff
        
#================================================================================================
from scipy.stats import pearsonr
def pearsonr2_TwoArComp(data1, data2, numbertype):
    """Returns r2 of two arrays for comparison"""
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]

    correlation, p_value = pearsonr(nonNA_data1, nonNA_data2)
    rsquare = correlation**2
    
    # calculate the difference
    Meandiff = np.nanmean(nonNA_data1-nonNA_data2)
    Meandiff_percent = Meandiff/np.nanmean(nonNA_data1)*100
    Mediandiff = np.nanmedian(nonNA_data1-nonNA_data2)
    
    # with specified format, print the results
    if numbertype=='f':
        print(r"r2 = %5.2f" % (rsquare))
        print(r"Mean Diff data1 - data2 = %5.2f; %5.0f percent" % (Meandiff,Meandiff_percent))
        print(r"Median of Diff data1 - data2 = %5.2f" % (Mediandiff))
    elif numbertype=='e':
        print(r"r2 = %5.2e" % (rsquare))
        print(r"Mean Diff data1 - data2 = %5.2e; %5.0f percent" % (Meandiff,Meandiff_percent))
        print(r"Median of Diff data1 - data2 = %5.2e" % (Mediandiff))

    return rsquare,Meandiff,Meandiff_percent,Mediandiff

def pearsonr_TwoArComp_returnr(data1, data2):
    """Returns r2 of two arrays for comparison"""
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]
    correlation, p_value = pearsonr(nonNA_data1, nonNA_data2)
    return correlation

#================================================================================================
from scipy.stats import spearmanr
def spearmanr2_TwoArComp(data1, data2, numbertype):
    """Returns r2 of two arrays for comparison"""
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]

    correlation, p_value = spearmanr(nonNA_data1, nonNA_data2)
    rsquare = correlation**2
    
    # calculate the difference
    Meandiff = np.nanmean(nonNA_data1-nonNA_data2)
    Meandiff_percent = Meandiff/np.nanmean(nonNA_data1)*100
    Mediandiff = np.nanmedian(nonNA_data1-nonNA_data2)
    
    # with specified format, print the results
    if numbertype=='f':
        print(r"r2 = %5.2f" % (rsquare))
        print(r"Mean Diff data1 - data2 = %5.2f; %5.0f percent" % (Meandiff,Meandiff_percent))
        print(r"Median of Diff data1 - data2 = %5.2f" % (Mediandiff))
    elif numbertype=='e':
        print(r"r2 = %5.2e" % (rsquare))
        print(r"Mean Diff data1 - data2 = %5.2e; %5.0f percent" % (Meandiff,Meandiff_percent))
        print(r"Median of Diff data1 - data2 = %5.2e" % (Mediandiff))

    return rsquare,Meandiff,Meandiff_percent,Mediandiff

def spearmanr_TwoArComp_returnr(data1, data2):
    """Returns r2 of two arrays for comparison"""
    # drop na values if appears as nan in either array
    nonNA_data1 = data1[(~np.isnan(data1))&~np.isnan(data2)]
    nonNA_data2 = data2[(~np.isnan(data1))&~np.isnan(data2)]
    correlation, p_value = spearmanr(nonNA_data1, nonNA_data2)
    return correlation

#================================================================================================
### regression: calculate the slope ###
import numpy as np
from scipy import stats

def regress2_slope(data1, data2, numbertype, method_type_2="ordinary least squares"):
    """
    Perform linear regression between two datasets (data1 and data2).

    Parameters:
        - data1: First dataset (numpy array or list)
        - data2: Second dataset (numpy array or list)
        - method_type_2: Regression method type ("ordinary least squares" or "reduced major axis")

    Returns:
        - slope: Slope of the regression line
        - intercept: Intercept of the regression line
    """
    # Convert data to numpy arrays
    data1 = np.array(data1)
    data2 = np.array(data2)

    # Check if method_type_2 is "ordinary least squares"
    if method_type_2.lower() == "ordinary least squares":
        slope, intercept, r_value, p_value, std_err = stats.linregress(data1, data2)
    elif method_type_2.lower() == "reduced major axis":
        # Calculate RMA regression
        # Assuming errors are approximately equal for both variables
        slope = np.std(data2) / np.std(data1) * np.sign(np.corrcoef(data1, data2)[0, 1])
        intercept = np.mean(data2) - slope * np.mean(data1)
    else:
        raise ValueError("Invalid method_type_2. Use 'ordinary least squares' or 'reduced major axis'.")

    print('slope: %5.2f' % slope)
    # with specified format, print the intercept
    if numbertype=='f':
        print('intercept: %5.2f' % intercept)
    elif numbertype=='e':
        print('intercept: %5.2e' % intercept)
    return slope, intercept

# # Example usage:
# data1 = [1, 2, 3, 4, 5]
# data2 = [2, 4, 5, 4, 5]

# # Calculate RMA regression slope
# rma_slope = regress2_slope(data1, data2, 'f', method_type_2="reduced major axis")
# # Print the RMA slope
# print("Reduced Major Axis (RMA) Regression Slope:", rma_slope)


#================================================================================================
### Compare model and observatopn using Reduced Major Axis Regression ###
import numpy as np

def model_vs_obs_RMA(Model_ar, Obs_ar, givenpval, numbertype):
    """
    Compare model vs observation using Reduced Major Axis (RMA) Regression.

    Parameters:
        - Model_ar (array-like): Array containing model data.
        - Obs_ar (array-like): Array containing observation data.

    Returns:
        - R2 (float): R-squared value for RMA Regression.
        - rma_slope (float): Slope of the RMA Regression line.

    Notes:
        - The R-squared value (R2) indicates the goodness of fit between the model and observation.
        - The slope (rma_slope) provides information about the relationship between the model and observation.
        - A slope < 1 suggests model underestimation, while a slope > 1 suggests model overestimation.
        - The function checks for NaN and inf values in both input arrays and removes them for analysis.
        - It also performs a Wilcoxon test to verify higher/lower values between model and observation.

    Example:
        >>> Model_ar = [2.0, 4.0, 5.0, 4.0, 5.0, np.nan, 8.0]
        >>> Obs_ar = [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0]
        >>> model_vs_obs_RMA(Obs_ar, Model_ar)
    """
    Datestr1 = "Model" 
    Datestr2 = "Observe"
    
    # Remove NaN and inf values in both arrays
    nonNA_Obs_ar = Obs_ar[(~np.isnan(Obs_ar)) & (~np.isnan(Model_ar))]
    nonNA_Model_ar = Model_ar[(~np.isnan(Obs_ar)) & (~np.isnan(Model_ar))]
    
    # Calculate R2
    R2 = ReducedMajorAxisRegression_R2_TwoArComp(nonNA_Model_ar, nonNA_Obs_ar)
    
    # Perform Wilcoxon test to verify higher/lower values
    # Model minus Observation (Datestr1-Datestr2)
    Meandiff,Meandiff_percent,Mediandiff = Oneside_Wilcoxon_test_result(Datestr1, nonNA_Model_ar, Datestr2, nonNA_Obs_ar, givenpval, numbertype)
    
    # Calculate the RMA regression slope (on the x-axis)
    rma_slope,rma_intercept = regress2_slope(nonNA_Obs_ar, nonNA_Model_ar, numbertype, method_type_2="reduced major axis")
    
    return R2, rma_slope, rma_intercept, Meandiff, Meandiff_percent, Mediandiff

# # example
# results = model_vs_obs_RMA(SurfObs_ar,CMAQ_ar, 0.05, 'f')

#================================================================================================
def OnlyReturnR2Slope_model_vs_obs_RMA(Model_ar, Obs_ar, givenpval):
    """
    Compare model vs observation using Reduced Major Axis (RMA) Regression.

    Parameters:
        - Model_ar (array-like): Array containing model data.
        - Obs_ar (array-like): Array containing observation data.

    Returns:
        - R2 (float): R-squared value for RMA Regression.
        - rma_slope (float): Slope of the RMA Regression line.

    Notes:
        - The R-squared value (R2) indicates the goodness of fit between the model and observation.
        - The slope (rma_slope) provides information about the relationship between the model and observation.
        - A slope < 1 suggests model underestimation, while a slope > 1 suggests model overestimation.
        - The function checks for NaN and inf values in both input arrays and removes them for analysis.

    Example:
        >>> Model_ar = [2.0, 4.0, 5.0, 4.0, 5.0, np.nan, 8.0]
        >>> Obs_ar = [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0]
        >>> model_vs_obs_RMA(Obs_ar, Model_ar)
    """
    Datestr1 = "Model" 
    Datestr2 = "Observe"
    
    # Remove NaN and inf values in both arrays
    nonNA_Obs_ar = Obs_ar[(~np.isnan(Obs_ar)) & (~np.isnan(Model_ar))]
    nonNA_Model_ar = Model_ar[(~np.isnan(Obs_ar)) & (~np.isnan(Model_ar))]
    
    # Calculate R2
    results = regress2(nonNA_Model_ar, nonNA_Obs_ar, _method_type_2="reduced major axis")
    R2 = np.square(results['r'])
    
    # Calculate the RMA regression slope (on the x-axis)
    rma_slope, rma_intercept, r_value, p_value, std_err = stats.linregress(nonNA_Obs_ar, nonNA_Model_ar)
    
    return R2, rma_slope, rma_intercept

#================================================================================================
def calculate_NRMSE(Model_ar, Obs_ar):
    """
    Calculate the Normalized Root Mean Square Error (NRMSE) between Model_ar and Obs_ar.
    NRMSE is a variation of RMSE that is normalized by the range of the observed values. 
    NRMSE values are unitless and range from 0 to infinity, with lower values indicating better model performance.

    Parameters:
    Model_ar (array-like): Array of model values.
    Obs_ar (array-like): Array of observed values.

    Returns:
    float: Normalized Root Mean Square Error (NRMSE).
    """
    # Calculate RMSE
    RMSE = np.sqrt(np.nanmean((Model_ar - Obs_ar) ** 2))
    
    # Calculate range of observed values
    obs_range = np.nanmax(Obs_ar) - np.nanmin(Obs_ar)
    
    # Calculate NRMSE
    NRMSE = RMSE / obs_range
    
    return NRMSE

#================================================================================================
# # Set the negative values to display lightblue
import pandas as pd
def color_negative(val):
    """
    Returns a CSS style string to set the background color based on the value.

    Parameters:
    val (numeric): The value to determine the background color.

    Returns:
    str: A CSS style string to set the background color.
    """
    if pd.isnull(val):  # Check for NaN values
        return ''
    try:
        color = 'lightblue' if float(val) < 0 else 'white'
    except ValueError:
        return ''  # Return empty string for non-numeric values
    return f'background-color: {color}'

# # styled_df = df.style.applymap(color_negative)


#================================================================================================
# def calculate_NMBE(Model_ar, Obs_ar):
#     """
#     Calculate the Normalized Mean Bias Error (NMBE).
    
#     Parameters:
#     Model_ar (array-like): Array or list of model predicted values.
#     Obs_ar (array-like): Array or list of observed values.
    
#     Returns:
#     float: Normalized Mean Bias Error (NMBE).
#     """
#     # Convert inputs to numpy arrays
#     Model_ar = np.array(Model_ar)
#     Obs_ar = np.array(Obs_ar)
    
#     # Calculate the bias error
#     bias_error = np.mean(Model_ar - Obs_ar)
    
#     # Calculate the mean of observed values
#     mean_observed = np.mean(Obs_ar)
    
#     # Calculate the NMBE
#     nmbe = bias_error / mean_observed
    
#     return nmbe

import numpy as np
def calculate_NMBE(Model_ar, Obs_ar):
    """
    Calculate the Normalized Mean Bias Error (NMBE).

    Parameters:
    Model_ar (array-like): Array or list of model predicted values.
    Obs_ar (array-like): Array or list of observed values.

    Returns:
    float: Normalized Mean Bias Error (NMBE) or None if unable to calculate.
    """
    # Convert inputs to numpy arrays
    Model_ar = np.array(Model_ar)
    Obs_ar = np.array(Obs_ar)

    # Mask out nan and inf in Model_ar and Obs_ar
    valid_indices = np.isfinite(Model_ar) & np.isfinite(Obs_ar)
    Model_ar = Model_ar[valid_indices]
    Obs_ar = Obs_ar[valid_indices]

    # Calculate the bias error
    bias_error = np.mean(Model_ar - Obs_ar)

    # Check if bias_error is not nan or inf
    if not np.isnan(bias_error) and not np.isinf(bias_error):
        # Calculate the mean of observed values
        mean_observed = np.mean(Obs_ar)

        # Check if mean_observed is not zero to avoid division by zero
        if mean_observed != 0:
            # Calculate the NMBE
            nmbe = bias_error / mean_observed
            return nmbe
        else:
            print("Error: Mean of observed values is zero, unable to calculate NMBE.")
    else:
        print("Error: Bias error contains nan or inf values, unable to calculate NMBE.")

# # Example usage:
# Model_ar = [11, 13, np.nan, 17, 19]
# Obs_ar = [10, 12, 15, 18, 20]

# nmbe_result = calculate_NMBE(Model_ar, Obs_ar)
# if nmbe_result is not None:
#     print("Normalized Mean Bias Error (NMBE):", nmbe_result)

#================================================================================================

import numpy as np

def calculate_MBE(Model_ar, Obs_ar):
    """
    Calculate the Mean Bias Error (MBE).

    Parameters:
    Model_ar (array-like): Array or list of model predicted values.
    Obs_ar (array-like): Array or list of observed values.

    Returns:
    float: Mean Bias Error (MBE) or None if unable to calculate.
    """
    # Convert inputs to numpy arrays
    Model_ar = np.array(Model_ar)
    Obs_ar = np.array(Obs_ar)

    # Mask out nan and inf in Model_ar and Obs_ar
    valid_indices = np.isfinite(Model_ar) & np.isfinite(Obs_ar)
    Model_ar = Model_ar[valid_indices]
    Obs_ar = Obs_ar[valid_indices]

    # Calculate the bias error
    bias_error = np.mean(Model_ar - Obs_ar)

    # Check if bias_error is not nan or inf
    if not np.isnan(bias_error) and not np.isinf(bias_error):
        return bias_error
    else:
        print("Error: Bias error contains nan or inf values, unable to calculate MBE.")
        return None
# # Example usage
# Model_ar = [1.2, 2.3, 3.1, 4.0, np.nan, 5.5]
# Obs_ar = [1.0, 2.1, 3.0, 4.1, 4.9, 5.4]
# MBE = calculate_MBE(Model_ar, Obs_ar)

import numpy as np

def calculate_relativeMBE(Model_ar, Obs_ar):
    """
    Calculate the Relative Mean Bias Error (rMBE) as:
        (mean(Model) - mean(Obs)) / mean(Obs) * 100

    Parameters:
    Model_ar (array-like): Array or list of model predicted values.
    Obs_ar (array-like): Array or list of observed values.

    Returns:
    float: Relative Mean Bias Error (%) or None if invalid.
    """
    import numpy as np

    # Convert inputs to numpy arrays
    Model_ar = np.array(Model_ar)
    Obs_ar = np.array(Obs_ar)

    # Mask invalid values
    valid_indices = np.isfinite(Model_ar) & np.isfinite(Obs_ar)
    Model_ar = Model_ar[valid_indices]
    Obs_ar = Obs_ar[valid_indices]

    if len(Model_ar) == 0 or len(Obs_ar) == 0:
        return None

    # Compute means
    mean_model = np.mean(Model_ar)
    mean_obs = np.mean(Obs_ar)

    # Relative MBE
    rel_bias = (mean_model - mean_obs) / mean_obs * 100

    if np.isfinite(rel_bias):
        return rel_bias
    else:
        return None



def calculate_RMSE(Model_ar, Obs_ar):
    """
    Calculate the Root Mean Square Error (RMSE) between Model_ar and Obs_ar.

    Parameters:
    Model_ar (array-like): Array of model values.
    Obs_ar (array-like): Array of observed values.

    Returns:
    float: Root Mean Square Error (RMSE).
    """
    # Convert inputs to numpy arrays
    Model_ar = np.array(Model_ar)
    Obs_ar = np.array(Obs_ar)

    # Mask out nan and inf in Model_ar and Obs_ar
    valid_indices = np.isfinite(Model_ar) & np.isfinite(Obs_ar)
    Model_ar = Model_ar[valid_indices]
    Obs_ar = Obs_ar[valid_indices]

    # Calculate RMSE
    RMSE = np.sqrt(np.mean((Model_ar - Obs_ar) ** 2))
    
    return RMSE

import numpy as np

def calculate_relativeRMSE(Model_ar, Obs_ar):
    """
    Calculate the Relative Root Mean Square Error (rRMSE, %) as:
        (RMSE / mean(Obs)) * 100

    Parameters:
    Model_ar (array-like): Array or list of model predicted values.
    Obs_ar (array-like): Array or list of observed values.

    Returns:
    float: Relative RMSE (%) or None if invalid.
    """
    # Convert inputs to numpy arrays
    Model_ar = np.array(Model_ar)
    Obs_ar = np.array(Obs_ar)

    # Mask out nan and inf
    valid_indices = np.isfinite(Model_ar) & np.isfinite(Obs_ar)
    Model_ar = Model_ar[valid_indices]
    Obs_ar = Obs_ar[valid_indices]

    if len(Model_ar) == 0 or len(Obs_ar) == 0:
        return None

    # RMSE
    RMSE = np.sqrt(np.mean((Model_ar - Obs_ar) ** 2))

    # Relative RMSE
    mean_obs = np.mean(Obs_ar)
    if mean_obs == 0:
        return None  # avoid division by zero
    rRMSE = RMSE / mean_obs * 100

    return rRMSE

# # Example usage
# Model_ar = [1.2, 2.3, 3.1, 4.0, np.nan, 5.5]
# Obs_ar = [1.0, 2.1, 3.0, 4.1, 4.9, 5.4]
# RMSE = calculate_RMSE(Model_ar, Obs_ar)

#================================================================================================


#================================================================================================
