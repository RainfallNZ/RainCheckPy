# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:41:10 2021

@author: Tim Kerr, Rainfall.NZ
"""

#Load modules
import numpy as np
import pandas as pd
import itertools
import math

def affinity( TestData, ReferenceData):
    'Compare the data between two sites to see how similar they are'
    'this uses an "affinity" index from Lewis et al. 2018, supplementary material'
    
    #join the two sets of data discarding periods not common to both
    result = TestData.join(ReferenceData, how='inner',lsuffix='_Test',rsuffix='_ref')
    result.columns =['Test','Reference']
    
    #Find the wet/dry values for each gauge
    TestWetDry = result.Test > 0
    ReferenceWetDry = result.Reference > 0
    
    #Combine the WetDry series together. Both wet = 2, both dry = 0
    BothWet = (TestWetDry * 1) + (ReferenceWetDry * 1)
    
    #Count the total number of both wets one wet or both drys
    CombinedTotals = BothWet.groupby(BothWet.values).count()
    #breakpoint()
    if ((0 in CombinedTotals) & (2 in CombinedTotals)):
        Affinity = CombinedTotals[0] / CombinedTotals.sum() + CombinedTotals[2] / CombinedTotals.sum()
    else:
        Affinity = 0
    return Affinity

def spearman( TestData, ReferenceData):
    "calculate the Spearman rank correlation coefficient between sites"
    
    #join the two sets of data discarding periods not common to both
    result = TestData.join(ReferenceData, how='inner',lsuffix='_Test',rsuffix='_ref')
    result.columns =['Test','Reference']
    
    CorrelationMatrix = result.corr(method="spearman")
    Spearman = CorrelationMatrix.Test['Reference']
        
    return Spearman
 
def neighborhoodDivergence( TestData: pd.DataFrame, ReferenceData: pd.DataFrame) -> pd.DataFrame:
    """Compares rainfall amounts to a another site
    
    Finds the ratio between the daily rainfall difference and the ninety-fifth percentile of
    the distribution of daily differences. This is analogous to the rain_outliers test
    but is based on comparison to an alternative site.
    This generates two values, the high divergence and the low divergence.
    High divergence is for when the Test value is higher than the reference value i.e. where the ratio of the max(0,Test - Reference) / 95th(max(0,Test - Reference)
    Low divergence is for when the Test value is lower than the reference value, i.e. ratio of the min(0,Test - Reference) / 5th(min(0,Test - Reference)
    
    Parameters
    ----------
    TestData : pd.DataFrame
        A time series of rainfall amounts for the site being tested.
    ReferenceData : pd.DataFrame
        A time series of rainfall amounts for the site to be compared with.

    Returns
    -------
    neighborhoodDivergence : pd.DataFrame
        A time series of 'LowOutlierData' and 'HighOutlierData'.
        High divergence is for when the Test value is higher than the reference value 
        i.e. where the ratio of the max(0,Test - Reference) / 95th(max(0,Test - Reference)
        Low divergence is for when the Test value is lower than the reference value, 
        i.e. ratio of the min(0,Test - Reference) / 5th(min(0,Test - Reference)

    """
    
    #Subset the TestData to just those observations of rain
    #TestDataNoRain = TestData.loc[TestData.Rainfall > 0,]

    #join the two sets of data discarding periods not common to both
    result = TestData.join(ReferenceData, how='inner',lsuffix='_Test',rsuffix='_ref')
    result.columns =['Test','Reference']
    
    #discard any dates with nan in either column
    result = result.dropna()
        
    #Create a series of absolute differences
    result['Differences'] = result.Test - result.Reference
     
    #Get the 95th percentile of the differences
    PosNinetyFifth = np.quantile(result.Differences[result.Differences > 0],0.95)
    NegFifth       = np.quantile(result.Differences[result.Differences < 0],0.05)

    #Find the ratio of each observation to the positive differences 95th percentile
    result['HighOutlierData'] = np.round(result.Differences / PosNinetyFifth,1)
    result.loc[result['HighOutlierData']<=0,'HighOutlierData'] = 0
    
    
    #Find the ratio of each observation to the negative differences 5th percentile
    result['LowOutlierData'] = np.round(result.Differences / NegFifth,1)
    result.loc[result['LowOutlierData']<=0,'LowOutlierData'] = 0
      
    #Join with the original data to include all the zero and nan observations
    neighborhoodDivergence = pd.merge(result,TestData,on='DateTime', how='right')[['LowOutlierData','HighOutlierData']]
    #neighborhoodDivergence = pd.merge(result,TestData,on='DateTime', how='right').fillna(0)[['LowOutlierData','HighOutlierData']]
    
       
    return neighborhoodDivergence

def DrySpellDivergence( TestData, ReferenceData):
    "find the ratio between the 15-day dry spell proportion difference and the ninety-fifth percentile of"
    "the distribution of the 15-day dry-spell proportion differences"
    
    #join the two sets of data 
    result = TestData.join(ReferenceData, how='outer',lsuffix='_Test',rsuffix='_ref')
    result.columns =['Test','Reference']
    
    #Restrict to when there is a time overlap
        #Get rid of the NaN's at the begining and end
    first_idx = max(TestData.first_valid_index(),ReferenceData.first_valid_index())
    last_idx = min(TestData.last_valid_index(),ReferenceData.last_valid_index())

    result = result.loc[first_idx:last_idx]
    
        #Add boolean columns for no-rain observations
    result['TestDry']=result.Test == 0
    result['ReferenceDry']=result.Reference == 0
    
    
    #Calculate the proportion of dry days in each 15 days for both gauges
    result['Test15dayDryCounts'] = result.TestDry.rolling(window='15d').sum()
    result['Reference15dayDryCounts'] = result.ReferenceDry.rolling(window='15d').sum()
    
    result['Test15dayObservationCounts'] = result.Test.rolling(window='15d').count()
    result['Reference15dayObservationCounts'] = result.Reference.rolling(window='15d').count()
    #Calculate a proportion difference series between the two gauges, but only when all 15 days were observed at both sites
    result['15DayDryProportionDifference'] = (result.Test15dayDryCounts / result.Test15dayObservationCounts) - (result.Reference15dayDryCounts / result.Reference15dayObservationCounts)
    
    #The following was replaced by the following following because it violates the python issue of setting a value on a copy of a slice!
    #result['15DayDryProportionDifference'][(result['Test15dayObservationCounts'] < 360)|(result['Reference15dayObservationCounts'] < 360)] = np.nan
    result.loc[(result['Test15dayObservationCounts'] < 360)|(result['Reference15dayObservationCounts'] < 360),'15DayDryProportionDifference'] = np.nan
    
    #Find the ninety fifth percentile of dry days portion positive differences (positive because we're only interested when the test is drier than the reference)
    if sum(result['15DayDryProportionDifference'] >= 0) > 0:
        DryProportionDiffereneNinetyFifth = np.quantile(result['15DayDryProportionDifference'][result['15DayDryProportionDifference'].notna() & (result['15DayDryProportionDifference'] >= 0)],0.95)
    else:
        DryProportionDiffereneNinetyFifth = 1
    #Calculate the ratio of the dry day proportion difference to the 95th percentile
    result['DryProportionOutlierIndex'] = np.round(result['15DayDryProportionDifference'] / DryProportionDiffereneNinetyFifth,1)
    result.loc[result['DryProportionOutlierIndex']<=0,'DryProportionOutlierIndex'] = 0
    
    #A big difference indicates suspect data
    DrySpellDivergence = result.DryProportionOutlierIndex
        
    return DrySpellDivergence

def TimeStepAllignment( TestData, ReferenceData):
    "This resamples the ReferenceData to match the observation times of the TestData"
    "this helps for comparison to irregularly sampled data (e.g. storage gauges"
    "or for manually recorded daily gauges that are read at non- 0:00 hours, e.g. at 8 or 9 am"
    ##For testing
    #TestData = LoadFromRainfallNZ_csv(RawDataFileName="../Data/PhD_data/Hooker Rd Bridge rainfall.csv")[0]
    #ReferenceData = LoadFromRainfallNZ_csv(RawDataFileName="../Data/PhD_data/Tasman Terminus Rainfall.csv")[0]
    
    #Merge the two timeseries. Add a third column which provides the aggregation index
    #join the two sets of data discarding periods not common to both
    result = TestData.join(ReferenceData, how='outer',lsuffix='_Test',rsuffix='_ref')
    result.columns =['Test','Reference']
        
    #Add a third column which provides the aggregation index
    result['aggregator'] = result.index.strftime('%Y-%m-%dT%H:%M%:%SZ')
    
    #replace all the aggregatior's without a Test observation with an nan
    result.aggregator[result.Test.isna()] = np.nan
    
    #interpolate the aggregator nan's with the following time value
    result['aggregator'] = result['aggregator'].fillna(method='bfill')
    
    #Perform a df.groupby of the reference data using the aggregation index
    g = result.groupby('aggregator')

    #Sum the groups, note the need to use aggregate and lambda instead of sum.
    #This is because sum(skipna=False) doesn't work.
    Reference_sums = g[['Reference']].aggregate(lambda x: sum(x))
    
    #Get rid of the NaN's at the begining and end
    first_idx = Reference_sums.first_valid_index()
    last_idx = Reference_sums.last_valid_index()

    Reference_sums = Reference_sums.loc[first_idx:last_idx]
    
    return Reference_sums

