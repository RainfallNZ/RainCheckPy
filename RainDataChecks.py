# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:41:10 2021

@author: Tim Kerr, Rainfall.NZ
"""

#Load modules
import numpy as np
import pandas as pd
import itertools
import numbers
import math
from scipy.signal import find_peaks
import pytz
import random
import scipy.stats as st

def rain_outliers( rain_data: pd.DataFrame ) -> pd.DataFrame:
    """Generates an outlier index time series
    
    Finds the ratio between each value to the ninety-ninth percentile of non-zero values.

    Parameters
    ----------
    rain_data : pd.DataFrame
        A time series of rainfall amounts to be tested.

    Returns
    -------
    rain_outliers : pd.DataFrame
        A time series of outlier indices. 

    """
    "Rainfall quality check for outliers"
    
    if len(rain_data.index) < 100:
        Output = pd.DataFrame(np.nan, columns=['Outlier'],index=rain_data.index)
    else:
        #Remove all the 0 observations
        NonZeroRainData = rain_data.values[rain_data.values > 0.2]
    
        #Find the ninetyninth percentile of the data
        NinetyNinth = np.quantile(NonZeroRainData,0.99)
    
        #Find the ratio of each observation to the 99th perentile
        OutlierData = np.round(rain_data.values / NinetyNinth,1)
        
        Output = pd.DataFrame(OutlierData, columns=['Outlier'],index=rain_data.index)
    
    return Output

def impossibles ( rain_data,minimum_precision=float('nan')):
    "rainfall quality check for impossible values"
    
    #Test for non-numbers, but don't include nan's
    NotANumber = ~np.array([isinstance(item,numbers.Number) for item in rain_data.values[:,0]])
    
    #Test for numeric values less than 0
    Sub_zeros = rain_data.apply(pd.to_numeric, errors='coerce') < 0
        
    #Test for numeric values that are not multiples of the minimum precision
    if not(math.isnan(minimum_precision)) and (minimum_precision > 0):
        False_precision = rain_data.apply(pd.to_numeric, errors='coerce') % minimum_precision != 0
    else: False_precision = NotANumber
    
    ImpossibleData = Sub_zeros.Rainfall.to_numpy() | NotANumber | False_precision
        
    Output = pd.DataFrame(ImpossibleData, columns=['Impossible'],index=rain_data.index)
    
    return Output

def DateTimeIssues (rain_data):
    "rainfall quality check for duplicate date times"
    
    DateTimeDuplicated = rain_data.index.duplicated(keep = False)
    
    Output = pd.DataFrame(DateTimeDuplicated, columns=['DuplicateDateTimes'],index=rain_data.index)
    
    return Output

def HighFrequencyTipping (rain_data):
    "rainfall quality check for unlikely rapid tipping"
    "from Blekinsop et al. (2017) lambda sub k statistic"
    "This is only appropriate for raw tip-based data"
    
    #Create a timeseries of inter-tip times in seconds
    InterTipTimes = rain_data.index.to_series().diff().astype('timedelta64[s]')
    
    #Test data
    #InterTipTimes = pd.Series([8,7,3,4,2,1,2,3,2,4,1,2,3,8,7,900,2,1,3,2,1,3,2,4,1,2,3,9,7,8])
    
    #Pre-allocate an array for the error series
    HighFrequencyTips = np.zeros(len(InterTipTimes),dtype=bool)
    
    #Calculate the lambda sub k statistic time series. This is a measure of the rate of change of inter-tip times
    LambdaSubK = np.log(InterTipTimes / InterTipTimes.shift(1)).abs()
    
    #Identify ocurrences of lambda sub k greater than 5
    RapidTipRateChanges = LambdaSubK > 5
 
    #create boolean series of sub threshold interTipTimes
    SubThresholdInterTipTimesBoolean = InterTipTimes < 5
      
    #For each Lanbdasubk below the threshold, get the indices of the 
    # immediately-following sequence of sub 5s intertiptimes. If there 
    # are 10 or more of them, mark the output for those time steps as suspect.
    RapidTipRateChangeIndices = np.where(RapidTipRateChanges)
    if len(RapidTipRateChangeIndices[0]) > 0:
        for index in RapidTipRateChangeIndices:
            #Find how many of the subsequent inter-tip times are less than 5 s
            SubThresholdTripTime = SubThresholdInterTipTimesBoolean[index[0]]
            NoOfSubThresholdTripTimes = 0
            while SubThresholdTripTime:
                NoOfSubThresholdTripTimes = NoOfSubThresholdTripTimes + 1
                SubThresholdTripTime = SubThresholdInterTipTimesBoolean[index[0] + NoOfSubThresholdTripTimes]
            HighFrequencyTips[index[0]:(index+NoOfSubThresholdTripTimes)[0]]= True 

    Output = pd.DataFrame(HighFrequencyTips, columns=['HighFrequencyTips'],index=rain_data.index)  
    return Output

def DrySpells (rain_data):
    "rainfall quality check for dry spells"
    "identify the length (in days) of a dry spell that a no-rain observation is within"
    "Alternative method using runlength encoding"
    
    #Create a one dimensional boolean panda series of dry/not dry
    DryObservations = pd.DataFrame((rain_data.values == 0),columns=["Dry"],index=rain_data.index)
    
    #Create a run-length-encoded version of the data
    RLE = [(k, sum(1 for i in g)) for k,g in itertools.groupby(DryObservations['Dry'])]
    
    
    RunLengthCodes = [a_tuple[0] for a_tuple in RLE]
    #Get the end index of each run by cusum the rle totals
    RunLengths = [a_tuple[1] for a_tuple in RLE]
    
    #Get the indices of the start and finish of each run
    RunLengthEndIndices = np.cumsum(RunLengths)-1
    RunLengthStartIndices = np.insert(RunLengthEndIndices[0:-1]+1,0,0,axis=0)
    
    #Get the subset that is for the dry runs
    DryRunLengthEndIndices = RunLengthEndIndices[RunLengthCodes]
    DryRunLengthStartIndices = RunLengthStartIndices[RunLengthCodes]
    
    #Get the start dates and end dates of the run lengths
    DryRunEndDateTimes = DryObservations.index[DryRunLengthEndIndices]
    DryRunStartDateTimes = DryObservations.index[DryRunLengthStartIndices]
    
    DryRunTimeLength = (DryRunEndDateTimes - DryRunStartDateTimes).days
    
    #Initialise a nan series ready to be populated with dry run day lengths
    DryObservations['DrySpellDayLengths'] = np.nan
    #DryObservations['RunLengthCodes'] = np.nan
    #breakpoint()
    #populate the new columns with the dry run day lengths 
    #Note mix of iloc and column name indexing, needed to avoid setting value of copy of a slice
    DryObservations.iloc[DryRunLengthEndIndices,DryObservations.columns.get_loc('DrySpellDayLengths')] = DryRunTimeLength
    #DryObservations['RunLengthCodes'].iloc[RunLengthEndIndices] = RunLengths
    #fill the gaps with next valid value
    DryObservations.DrySpellDayLengths.fillna(method='backfill',inplace=True)
    #DryObservations.RunLengthCodes.fillna(method='backfill',inplace=True)
    #Make all non-zero rainfalls have run lengths of 0
    DryObservations.loc[~DryObservations.Dry,'DrySpellDayLengths'] = 0
    
    Output = DryObservations[['DrySpellDayLengths']]
    
    return Output


def RepeatedValues (rain_data):
    "rainfall quality check for unlikely repeating values"
    "identify the length (in consecutive time units) that a value is repeated"
    "this check should not be applied to tip data"
    
    #Create a one dimensional boolean panda series of wet/not wet
    WetObservations = pd.DataFrame(((rain_data.values > 0) * rain_data.values),columns=["Wet"],index=rain_data.index)
    
    #Create a run-length-encoded version of the data
    RLE = [(k, sum(1 for i in g)) for k,g in itertools.groupby(rain_data.iloc[:,0])]
    
    RunLengthCodes = np.array([a_tuple[0] for a_tuple in RLE])
    #Get the end index of each run by cusum the rle totals
    RunLengths = np.array([a_tuple[1] for a_tuple in RLE])
    
    #Get the indices of the finish of each run
    RunLengthEndIndices = np.cumsum(RunLengths)-1
    
    #Get the subset that is for the wet runs
    WetRunLengthEndIndices = RunLengthEndIndices[RunLengthCodes>0]
    WetRunLengths = RunLengths[RunLengthCodes>0]
    #WetRunLengthStartIndices = RunLengthStartIndices[RunLengthCodes]
    
    #Initialise a nan series ready to be populated with dry run day lengths
    WetObservations['RepeatedValues'] = np.nan
    #DryObservations['RunLengthCodes'] = np.nan
    
    #populate the new columns with the wet run lengths 
    WetObservations.iloc[WetRunLengthEndIndices,WetObservations.columns.get_loc('RepeatedValues')] = WetRunLengths
    #DryObservations['RunLengthCodes'].iloc[RunLengthEndIndices] = RunLengths
    #fill the gaps with next valid value
    WetObservations.RepeatedValues.fillna(method='backfill',inplace=True)
    #DryObservations.RunLengthCodes.fillna(method='backfill',inplace=True)
    #Make all non-zero rainfalls have run lengths of 0
    WetObservations.loc[~(WetObservations.Wet > 0),'RepeatedValues'] = 0
    
    Output = WetObservations[['RepeatedValues']]
    
    return Output

def Homogeneity (rain_data):
    """Applies the Pettitt non-parameteric test to annual series to determine if there are major inhomogeneities in the data
       If there is, the test is repeated on the most recent side of the inhomogeneity to test if there is another.
       The most recent section that is homogeneous is retained and the remainder flagged.
       This uses the pyHomogeneity package https://github.com/mmhs013/pyHomogeneity
    """
    import pyhomogeneity as hg

    if len(rain_data.index) < 100:
        Output = pd.DataFrame(np.nan, columns=['Homogeneous'],index=rain_data.index)
    else:
        #Initialise Homogeneous timeseries assuming everything is OK 
        Homogeneous = pd.DataFrame(True, columns=['Homogeneous','ChangePoint'],index=rain_data.index)
        
        #Create an annual series from the data if it exists, using years with at least 96 % of a year (i.e. 11 and a half months)
        DataStepLengthInHours = (rain_data.index[1] - rain_data.index[0]).total_seconds()//3600
        AnnualData = rain_data.resample("1y").sum(min_count = int(0.96 * 365 * 24 / DataStepLengthInHours ))
        #Check for homogeneity if there are more than 3 years of data
        if AnnualData.count().any() > 3:   
            #Apply the Pettitt test to the annual totals
            result = hg.pettitt_test(AnnualData)
            #While there is inhomogeneity, Keep checking the series after it
            MoreInhomogeneity = result.h
            while MoreInhomogeneity:
                Homogeneous[Homogeneous.index < pd.to_datetime(result.cp)]=False
                if len(AnnualData[AnnualData.index > pd.to_datetime(result.cp)]) > 3:
                    result=hg.pettitt_test(AnnualData[AnnualData.index > pd.to_datetime(result.cp)])
                    MoreInhomogeneity = result.h
                else: MoreInhomogeneity = False
    
        Output = Homogeneous
    return Output

def SubFreezingRain (rain_data, temperature_data):
    """"rainfall quality check for observations during freezing temperatures
    identify the observations when the maximum temperature was less than zero degrees C
    """

    #merge the rain and temperature data together
    RainAndTemperature = pd.merge(left = rain_data,right = temperature_data, left_index=True,right_index=True, how = 'left')

    #find when it is raining and the temperature is less than zero
    RainAndTemperature['FreezingRain'] = (RainAndTemperature.Rainfall > 0) & (RainAndTemperature.TMax < 0)  
    
    Output = RainAndTemperature['FreezingRain']
    
    return Output

def RelatedFlowEvents (rain_data, Daily_streamflow_data):
    """"rainfall quality check for observations compared to flow events
    for each time step allocate the relative magnitude of a peak flow event ocurring on the same day or the day after
    but only if rain events are associated with flow events
    used with daily streamflow and hourly rainfall, possibly daily rainfall, but it hasn't been tested yet.'
    """
    #For testing, get some rain data and some stream flow data
    #from LoadDataFunctions import LoadFrom_ClimateDataBase_netCDF
    #rain_data = LoadFrom_ClimateDataBase_netCDF(AgentNumber = 17610)['Hourly'] #This is Snowdon Raws site near the head of the Selwyn
    #rain_data = LoadFrom_ClimateDataBase_netCDF(AgentNumber = 41489)['Hourly'] #This is Arthurs Pass
    #Check the timezone, and if not set, set it to NZST
    #if rain_data.index.tzinfo is None or rain_data.index.tzinfo.utcoffset(rain_data.index) is None: rain_data.index=  rain_data.index.tz_localize(pytz.timezone("Etc/GMT-12"))
    
    #from LoadDataFunctions import LoadFromTethysDownloads_csv
    #Daily_streamflow_data = LoadFromTethysDownloads_csv()['Daily'] #This is Selwyn at Whitecliffs
    #Daily_streamflow_data = LoadFromTethysDownloads_csv(RawDataFileName ='../../Data/tethysDownloads/Waimakariri River at Otarama.csv')['Daily']
    
    #For testing, restrict to last 3 years
    #Daily_streamflow_data = Daily_streamflow_data.iloc[-1095:,]

    #Only apply test if there is more than two days of data
    if  max(rain_data.index)- min(rain_data.index) < pd.Timedelta('2 days'):
        RainAndFlow = rain_data.copy()
        RainAndFlow['Peak_prominence'] = np.nan
    else: 
        
        #Find flow peaks, where a peak is higher than the inter-peak low by at least 20 % of the mean flow.
        #This definition should identify most peaks without getting the tiny variations
        #the function returns (among other things) the prominence of each peak, which is the 
        #vertical difference between the peak and the lowest point within 'wlen' of the peak, 
        #or to the next peak that is higher than the current peak, if that is less than 'wlen'.
        
        peaks = find_peaks(Daily_streamflow_data['Streamflow'], height = 0, prominence=Daily_streamflow_data.mean().item() * 0.1,wlen = 3)
        DaysWithPeaks = Daily_streamflow_data.index[peaks[0]]
        #Get the 95th percentile of the peak prominences
        NinetyFifthPP = np.quantile(peaks[1]['prominences'],0.95)
        
        #Find the ratio of each peak's prominence to the 95th perentile
        RelativeProminence = np.round(peaks[1]['prominences'] / NinetyFifthPP,3)
        
        #Create a data frame with the same index as the flow data, but with the peak data in it on the days with peaks
        FlowPeakSeries = pd.DataFrame(index = Daily_streamflow_data.index, columns = ['Peak_prominence'])
        FlowPeakSeries.loc[DaysWithPeaks,'Peak_prominence'] = RelativeProminence
        
        #Strip timezone from FlowPeakSeries to make it compatible with the rain and temperature data
        if FlowPeakSeries.index.tzinfo is not None: 
            FlowPeakSeries.index=  FlowPeakSeries.index.tz_convert(pytz.timezone("Etc/GMT-12"))
            FlowPeakSeries.index=  FlowPeakSeries.index.tz_localize(None)
        
        #merge the rain and flow peak data together
        RainAndFlow = pd.merge(left = rain_data,right = FlowPeakSeries['Peak_prominence'], left_index=True,right_index=True, how = 'left')
        
        
        #fill flow peak so that the whole day has the flow peak relative prominence value, and backfill for 24 hours to associate with possible rain, accounting for time-to-concentration
        FillLength = int(24 // ((RainAndFlow.index[1] - RainAndFlow.index[0]).total_seconds()//3600))
        RainAndFlow.loc[:,'Peak_prominence'] = RainAndFlow.loc[:,'Peak_prominence'].fillna(method="pad",limit=FillLength-1).fillna(method='bfill', limit=FillLength).fillna(0)
        
        #Determine if high rainfall events are associated with peak flow events
        #GetRainEvents > 99th percentile
        #Get peak prominence for those events
        NonZeroRainData = rain_data.values[rain_data.values > 0]
    
        #Find the ninetyninth percentile of the data
        NinetyNinth = np.quantile(NonZeroRainData,0.99)
        
        #Find the rain events that are greater than the 99th
        HighRainHours = RainAndFlow['Rainfall'] > NinetyNinth
        
        #Need to get each rain events maximum hourly rainfall, where events are spearated by at least 12 hours.
        #Get the time difference between each high rainfall event
        TimeDifferenceSeriesOnhighRainEvents = HighRainHours[HighRainHours].index.to_series().diff().to_frame()
        
        #Round the differences to the nearest 6 hours
        TimeDifferenceSeriesOnhighRainEvents['RoundedDateTime'] = TimeDifferenceSeriesOnhighRainEvents['DateTime'].dt.round('12H')
        
        #Add the rain and the flow peak prominence values back on
        TimeDifferenceSeriesOnhighRainEvents['Rainfall']        = RainAndFlow.loc[HighRainHours,'Rainfall']
        TimeDifferenceSeriesOnhighRainEvents['Peak_prominence'] = RainAndFlow.loc[HighRainHours,'Peak_prominence']
        
        #Apply magic to group by RoundedTDateTime and get maximum rainfall and flow peak prominence for each event
        RainEventPeakProminence = TimeDifferenceSeriesOnhighRainEvents.groupby((TimeDifferenceSeriesOnhighRainEvents['RoundedDateTime'] != TimeDifferenceSeriesOnhighRainEvents['RoundedDateTime'].shift()).cumsum(), as_index=False).agg(
            {'DateTime': 'first', 'Rainfall': 'max', 'Peak_prominence': 'max'})
            
        #Get peak prominence for random hours
        NumberOfSamples = min(len(RainAndFlow.index),10000)
        RandomPeakProminence = RainAndFlow.loc[RainAndFlow.index[random.sample(range(0,len(RainAndFlow.index)),NumberOfSamples)],'Peak_prominence']
        
        #Likelihood of flow peak event from random hours
        RandomPeakLilelihood = np.count_nonzero(RandomPeakProminence)/NumberOfSamples
               
        #Test whether the likelihood of flow events during high rainfall events is different from the likelihood flow events during random hours.
        #Use a binomial test.
        #The null hypothesis is that the likelihoods match. So if the result is less than 0.01, this means there is a small likelihood that they are the same, so they can be considered different.
        #If the result is > 0.01 then they're considered the same, and so flow events are not going to be helpful in confirming high rainfall events.
        ProbabilityThatFlowEventsDuringHighRainEventsMatchesRandom = st.binom_test(x=np.count_nonzero(RainEventPeakProminence['Peak_prominence']),
                                                                                   n=RainEventPeakProminence['Peak_prominence'].count(),p=RandomPeakLilelihood)
        
        #If the flow peaks are not statistically related to rain events, reset the Peak_provenance to NA
        if ProbabilityThatFlowEventsDuringHighRainEventsMatchesRandom > 0.01:
           RainAndFlow['Peak_prominence'] = np.nan
    

    
    Output = RainAndFlow['Peak_prominence']
    
    return Output




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

    #Get the 95th percentile of the differences, but only if there are some positive differences
    if sum(result.Differences > 0) > 0:
        PosNinetyFifth = np.quantile(result.Differences[result.Differences > 0],0.95)
        #Find the ratio of each observation to the positive differences 95th percentile
        result['HighOutlierData'] = np.round(result.Differences / PosNinetyFifth,1)
        result.loc[result['HighOutlierData']<=0,'HighOutlierData'] = 0
    else:
       result['HighOutlierData'] = 0
     

    #Get the 5th percentile of differences, but only if there are some negative differences
    if sum(result.Differences < 0) > 0:
        NegFifth       = np.quantile(result.Differences[result.Differences < 0],0.05)
        #Find the ratio of each observation to the negative differences 5th percentile
        result['LowOutlierData'] = np.round(result.Differences / NegFifth,1)
        result.loc[result['LowOutlierData']<=0,'LowOutlierData'] = 0
    else:
       result['LowOutlierData'] = 0 
      
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
