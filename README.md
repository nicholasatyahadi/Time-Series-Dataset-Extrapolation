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
* Correlation between values, especially in multivariate dataset.
* Standard deviation or Margin of Error (MOE) of the dataset.

In this example, I will show you how to work with univariate (with and without MOE) and multivariate time series dataset. 
The dataset that I use in this example will be linked below.

### Univariate Dataset with MOE
### Univariate Dataset without MOE
### Multivariate Dataset
