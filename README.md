# Time-Series-Dataset-Extrapolation
## Overview
Sometimes we receive data samples with small amount of data. This small amount might be affected by the periodcity of the observation.
For example, daily observation will result in a larger set of data than annual observation.

Say that you need to create a model based on annual observations but the amount of observations is just too small to create an appropiate model.
One way to do it is to extrapolate the dataset by keeping the original properties of the dataset as best as we can.
With this in mind I've created one of many methods you can extrapolate your data in such manner.

## How to Extrapolate
Here are a few things to keep in mind while creating an extrapolation method for the dataset:
* Data characteristics: The trend or distribution of the dataset.
* Data cleaning: Making sure the dataset contains the appropiate values
* Correlation between values, especially in multivariate dataset.
* Standard deviation or Margin of Error (MOE) of the dataset.

In this example, I will show you how to work with univariate (with and without MOE) and multivariate time series dataset. 
The dataset that I use in this example will be linked below. Also, for this example I'm going to use Uniform(0,1) distribution as the distribution for the standardized values.
(Yes you can use other distributions if you want to)

### Univariate Dataset with MOE
We know that the MOE is the error amount of the random sampling of a survey. So knowing the amount of the MOE gives an advantage in understanding the dataset better.
First we're going to load the dataset, which is the US Annual Mean Income, gathered from United Census Bureau.
### Univariate Dataset without MOE
### Multivariate Dataset
