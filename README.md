# Time-Series-Dataset-Extrapolation
## Overview
Sometimes we receive data samples with small amount of data. This small amount might be affected by the periodcity of the observation. For example, daily observation will result in a larger set of data than annual observation.

Say that you need to create a model based on annual observations but the amount of observations is just too small to create an appropiate model. One way to do it is to extrapolate the dataset by keeping the original properties of the dataset as best as we can. With this in mind I've created one of many methods you can extrapolate your data in such manner.

## How to Extrapolate
Here are a few things to keep in mind while creating an extrapolation method for the dataset:
* Data characteristics: The trend or distribution of the dataset.
* Data cleaning: Making sure the dataset contains the appropiate values
* Correlation between values, especially in multivariate dataset.
* Standard deviation or Margin of Error (MOE) of the dataset.

In this example, I will show you how to work with univariate (with and without MOE) and multivariate time series dataset. 
The dataset that I use in this example will be linked below. Also, for this example I'm going to use Uniform(0,1) distribution as the distribution for the standardized values (yes you can use other distributions if you want to) and set the seed to 23.

The standardization equation:

<img width="178" alt="stand" src="https://user-images.githubusercontent.com/53423050/91012503-69603780-e610-11ea-936b-4eb2982eb7ea.png">

### Univariate Dataset with MOE
We know that the MOE is the error amount of the random sampling of a survey. So knowing the amount of the MOE gives an advantage in understanding the dataset better.
First we're going to load the dataset, which is the US Annual Mean Income, gathered from United States Census Bureau Table S1901: https://data.census.gov/

<img width="265" alt="income" src="https://user-images.githubusercontent.com/53423050/91012412-446bc480-e610-11ea-8719-bd46d5a5a09b.png">

I've already added the Min and Max columns which are derived by adding and substracting the mean with the MOE, these will serve as our minimum and maximum values to be used in standardization equation.

Say it's best to have monthly values from the dataset. Hence, we're going to generate 12 random uniform values to act as our standardized values for each year. These randomly generated values will be then converted back to the real income values. The randomization process will be adjusted according to the minimum and maximum values each available year. Here are the before and after plot of the dataset:

Original annual dataset:

<img width="271" alt="tsincomereal" src="https://user-images.githubusercontent.com/53423050/91012536-74b36300-e610-11ea-97dd-d6592eb67076.png">

Here's the extrapolated one:

<img width="262" alt="tsincomeex" src="https://user-images.githubusercontent.com/53423050/91012553-7c730780-e610-11ea-8026-0a71938f9477.png">

As we can see, the increasing trend of the income is still existent. Even so, the extrapolated values are the noisier version of the original values, which could be used to model a complex time series model. To add, here's the percentage difference of statistical properties between the two to show you how accurate it is:

<img width="167" alt="sumdiffincomoe" src="https://user-images.githubusercontent.com/53423050/91012580-84cb4280-e610-11ea-8c18-ae7f20e3400f.png">

Pretty close isn't it? :)

### Univariate Dataset without MOE
For this example we're going to use the Monthly Sunspots Dataset which can be downloaded here:
https://machinelearningmastery.com/time-series-datasets-for-machine-learning/

Also, since the original dataset is already large, we're only going to use the last 36 monthly observations as our example. Here is the plot of the data:

<img width="262" alt="tssunspotsreak" src="https://user-images.githubusercontent.com/53423050/91016428-f0b0a980-e616-11ea-82d8-30e1a31ea2c5.png">

Say that we're required to generate another set of monthly observations based on the extracted dataset.
### Multivariate Dataset
