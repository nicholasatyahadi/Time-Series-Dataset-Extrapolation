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

In this example, I will show you how to work with univariate (with and without MOE) and multivariate time series dataset. The dataset that I use in this example will be linked below. Also, for this example I'm going to use Uniform(0,1) distribution as the distribution for the standardized values (yes you can use other distributions if you want to) and set the seed to 23.

Here is the standardization equation:

<img width="300"  alt="stand" src="https://user-images.githubusercontent.com/53423050/91012503-69603780-e610-11ea-936b-4eb2982eb7ea.png">

To measure the accuracy of the extrapolation method, we're going to compare the statistical properties, which are minimum, maximum, average, then first, second, and third quantiles, of both real and extrapolated datasets. The values will be compared using the percentage difference between the values of statistical properties by this method:

<img width="275" alt="statprop" src="https://user-images.githubusercontent.com/53423050/91028037-18a80900-e627-11ea-9634-cc41e793565a.png">

### Univariate Dataset with MOE
We know that the MOE is the error amount of the random sampling of a survey. So knowing the amount of the MOE gives an advantage in understanding the dataset better.
First we're going to load the dataset, which is the US Annual Mean Income, gathered from United States Census Bureau Table S1901: https://data.census.gov/

<img width="350" alt="income" src="https://user-images.githubusercontent.com/53423050/91012412-446bc480-e610-11ea-8719-bd46d5a5a09b.png">

I've already added the Min and Max columns which are derived by adding and substracting the mean with the MOE, these will serve as our minimum and maximum values to be used in standardization equation.

Say it's best to have monthly values from the dataset. Hence, we're going to generate 12 random uniform values to act as our standardized values for each year. These randomly generated values will be then converted back to the real income values. The randomization process will be adjusted according to the minimum and maximum values each available year. Here are the before and after plot of the dataset:

Original annual dataset:

<img width="350" alt="tsincomereal" src="https://user-images.githubusercontent.com/53423050/91012536-74b36300-e610-11ea-97dd-d6592eb67076.png">

Here's the extrapolated one:

<img width="350" alt="tsincomeex" src="https://user-images.githubusercontent.com/53423050/91012553-7c730780-e610-11ea-8026-0a71938f9477.png">

As we can see, the increasing trend of the income is still existent. Even so, the extrapolated values are the noisier version of the original values, which could be used to model a complex time series model. To add, here's the percentage difference of statistical properties between the two to show you how accurate it is:

<img width="330" alt="sumdiffincomoe" src="https://user-images.githubusercontent.com/53423050/91012580-84cb4280-e610-11ea-8c18-ae7f20e3400f.png">

Pretty close isn't it? :)

### Univariate Dataset without MOE
For this example we're going to use the Monthly Sunspots Dataset which can be downloaded here:
https://machinelearningmastery.com/time-series-datasets-for-machine-learning/

Also, since the original dataset is already large, we're only going to use the last 36 monthly observations as our example. Here is the plot of the data:

<img width="350" alt="tssunspotsreak" src="https://user-images.githubusercontent.com/53423050/91024350-07a8c900-e622-11ea-821f-c8303f5d4e3e.png">

Say that we're required to generate daily observations based on the extracted dataset. Hence, we're going to need the MOE for the dataset, especially the MOE of the monthly values. One way to get the MOE is to fit the existing values to a model then retrive the parameters and the standard errors. Since we know that this is a monthly dataset, we can use linear seasonal regression as the model to be fitted. Moreover, by using regression we must assume that the dataset follows a homoscedasticity[1] assumption. Here's the result of the fit:

<img width="400" alt="seasonressun" src="https://user-images.githubusercontent.com/53423050/91023948-66217780-e621-11ea-853f-83b2469b5a4c.png">

Which is already a great fit. Now we're only going to use the monthly parameters (estimate) and standard errors to create the 95% confidence interval for each month. The confidence interval will be used to randomize values with the same method as explained in the last section.

Here's the plot for the daily extrapolated value:

<img width="550" alt="dailysunspots" src="https://user-images.githubusercontent.com/53423050/91026464-f8774a80-e624-11ea-8f3f-ed261f18507e.png">

Yes, it's super noisy. But let's look at the percentage difference of the statistical properties of each datasets

<img width="300" alt="sumdiffsun" src="https://user-images.githubusercontent.com/53423050/91026652-412f0380-e625-11ea-9fdd-02772f92fd10.png">

This time, it's not as impressive as the income example. Why? First, we can definitely see that the extrapolated values doesn't carry the decreasing trend in the values and is stationary[2]. Morever, since the extrapolated values are daily observations and the values are stationary, the standard deviation of the observations decreases, causing the percentage difference of -20%. These findings told us immediately that linear seasonal regression may not be the perfect fit for the dataset.

### Multivariate Dataset
For our last example, we're going to use the New Delhi Climate training dataset, which can be acquired here:

https://www.kaggle.com/sumanthvrao/daily-climate-time-series-data

Now here's the twist, since this is a multivariate dataset, we can't simply immediately find a model that would suffice the dataset's characteristics. Let's start observing step by step:

1. Data Validation

The dataset contains four climate properties: Temperature, Humidity, Wind Speed, and Air Pressure. Let's take a look at the scatterplot between those fields

<img width="463" alt="corrrawdaily" src="https://user-images.githubusercontent.com/53423050/91029108-8ef93b00-e628-11ea-8e80-d3d088ce28c3.png">

Besides the apparent relation between temperature and humidity, there's something off about the pressure values. Most of the values have similar values (as seen in the thick blue cluster) and around 1000, but there are definitely outliers as shown in the small dots which are away from the cluster. Moreover, if we summarize the values by month and observe the minimum and maximum pressures, there are definitely a few outliers.

<img width="360" alt="weirdmin" src="https://user-images.githubusercontent.com/53423050/91031045-3b87ec80-e62a-11ea-8ff4-dc84049ca567.png">
<img width="360" alt="weirdmax" src="https://user-images.githubusercontent.com/53423050/91031050-3cb91980-e62a-11ea-95e8-af67c39e6466.png">

I am no climate expert, so why don't we ask our smart best friend, jack-of-all-trades, Google. Based on the facts that I've found in https://iridl.ldeo.columbia.edu/dochelp/QA/Basic/atmos_press.html the average air pressure is around 1013.25 millibars. Hence, with that in mind along with the values in the dataset, we can determine that air pressure between 990 and 1024 is normal, the rest is outliers. The outliers will be adjusted using the distribution of values according to the month.

Here are the average and standard deviation summary of each month after cleaning:

<img width="700" alt="climsumclean" src="https://user-images.githubusercontent.com/53423050/91031743-37a89a00-e62b-11ea-95db-12efafa6a0eb.png">

Never been better! Now here's the new scatterplots:

<img width="452" alt="cleancorr" src="https://user-images.githubusercontent.com/53423050/91031858-560e9580-e62b-11ea-86c2-58413dad9aae.png">


2. Data Characteristics

Now we're off to finding characteristics between the plot. BUT HERE'S A TWIST!

You didn't actually receive the daily values, instead, you were given the monthly average values and asked if it's possible to create daily values from it? (Try to guess the answer and see the results below, keep reading!)

Let's take a look at the correlation of the monthly average dataset.

Scatterplots:

<img width="451" alt="avgcorr" src="https://user-images.githubusercontent.com/53423050/91032578-59565100-e62c-11ea-91ac-bfc68a450794.png">

Correlation plot:

<img width="462" alt="corrplot" src="https://user-images.githubusercontent.com/53423050/91032851-beaa4200-e62c-11ea-9a44-37d80c65ae62.png">

Here's what we can see from the plot: There are correlations between the climate properties ranging from moderate to high correlation. (sigh..)

But let's take a closer look in the scatterplots. The temperature and pressure relation definitely shows a negative correlation with quadratic pattern, whereas the other relations seems to only show negative/positive correlation with unclear pattern. With this in mind, why don't we create a linear and quadratic regression between the fields and decide which pair of fields are actually effectively correlating to each other using the R-Squared values? But before doing so, we need to standardized the values in each field since the values from each field is not the same.

Here are the model fitting results along with the result of fitting each field to linear seasonal regression model:

<img width="400" alt="regdf" src="https://user-images.githubusercontent.com/53423050/91033913-59575080-e62e-11ea-84c1-2aa84271bd66.png">

Linearr2 and quadr2 are the r-squared values from the linear and quadratic regression fit respectively. While the adjlinearr2 and adjquadr2 are the adjusted r-squred values from the linear and quadratic regression fir respectively. The quadr2 and adjquadr2 for the seasonal regression is set to 0 since we didn't do any quadratic regression fit to it.

By eliminating the results from the seasonal trend, we can see from the table the highest r-squared values with pressure as the dependent variable and temperature as the independent variable. Hence, when we simulate the temperature we need to simulate the pressure values first then use the quadratic model to acquire the temperature values. Meanwhile, the other regression results don't seem quite significant despite the high correlations we found earlier. Therefore, we will assume that the only relation that exists are pressure affects the temperature and the rest are independent to each other. To sum it all up, here's what are we going to do:

<img width="572" alt="flow" src="https://user-images.githubusercontent.com/53423050/91185903-926bef80-e718-11ea-8041-708d6add9fb6.png">

Note that to extrapolate the temperature values, we are using the standard error geneerated from the seasonal model for temperature. With that in mind and the workflow has been done, here are plots comparing the real vs extrapolated values (the real values shown with the black graph):

<img width="500" alt="avgtemp" src="https://user-images.githubusercontent.com/53423050/91062278-b7982980-e656-11ea-9f2a-1f160d4aa6d2.png">

<img width="500" alt="avghumid" src="https://user-images.githubusercontent.com/53423050/91062303-bbc44700-e656-11ea-9139-594426275a9f.png">

<img width="500" alt="avgspeed" src="https://user-images.githubusercontent.com/53423050/91062273-b535cf80-e656-11ea-92e1-77ec61419544.png">

<img width="500" alt="avgpress" src="https://user-images.githubusercontent.com/53423050/91062308-bc5cdd80-e656-11ea-8958-0728e3336521.png">

To compare the values, I'm going to introduce you to another way to compare, which is using the average of percentage difference by time points. The percentage difference is still calculated the same way as the statistical properties. Here is the result

![image](https://user-images.githubusercontent.com/53423050/91037511-53b03980-e633-11ea-8c89-d2af87a1eadc.png)

After looking at the comparison plot and the average percentage differences, we can consider the extrapolation method is quite appropiate. Unfortunately, the real question has not been answered! :)

We still need to convert it back to daily values then comparing it. Since we already has the standard errors gathered from each seasonal model, we can use these values to estimate the margin of errors for the daily values for each field. But how so? Let me get a little bit mathematical for a bit.

Let Xi be the daily values and Yj be the monthly average. We also know that margin of error is defined by

<img width="107" alt="moe" src="https://user-images.githubusercontent.com/53423050/91038095-69722e80-e634-11ea-9287-0c3b9fd2660b.png">

in which Z is the Normal distribution deviate according to alpha and sigma/sqrt(n) is the standard error. Since we know the standard error, we can find the standard deviation by multiplying the standard error by sqrt(n), in which n is the number of monthly observations. Do note that the variance is for the standardized values so we need to convert the variance by these following set of equations:

<img width="350" alt="varyj" src="https://user-images.githubusercontent.com/53423050/91053630-4bfd8e80-e64d-11ea-94d4-799ada2cccee.png">

Then we will have the estimation of the standard error from the daily values, using these set of equations:

<img width="350" alt="varxi" src="https://user-images.githubusercontent.com/53423050/91053633-4d2ebb80-e64d-11ea-89cc-2bdd59ff7b80.png">

in which nj is the number of days in month j.

Now with that been said, we can extrapolate the monthly values into the daily values. Let's see the comparison plots:

<img width="500" alt="tempdailiy" src="https://user-images.githubusercontent.com/53423050/91062297-bb2bb080-e656-11ea-958d-d03555e75153.png">

<img width="500" alt="dailyhumid" src="https://user-images.githubusercontent.com/53423050/91062283-b7982980-e656-11ea-9764-17c346f49385.png">

<img width="500" alt="dailywind" src="https://user-images.githubusercontent.com/53423050/91062293-b9fa8380-e656-11ea-9398-71a9bb5a0ca0.png">

<img width="500" alt="dailypress" src="https://user-images.githubusercontent.com/53423050/91062285-b830c000-e656-11ea-9c19-d1e1373df249.png">

I don't think we need to deep dive to calculate the average percentage difference for each field. It's apparent that the extrapolation method didn't replicate the actual daily values, here's why:

* After the conversion, the dataset becomes too noisy and the extrapolation method only able to replicate for a short-term.
* Since the average values doesn't have a "jump" pattern especially in wind speed, the extrapolated daily values won't be able to catch that.

So, the answer is definitely a no.

### Limitations of the Method
The extrapolated values from a dataset can only replicate the characteristics of the dataset used in the extrapolation process. Moreover, by assuming a homocedasticity condition in the calculations may be the cause of the inability to detect "jump"(s). Therefore, an appropiate time-series model to fit the values are important to acknowledge to get better results. But let's not forget, the extrapolated values won't be and never be a perfect representation of the population. Nevertheless, it's still useful to estimate and understand what the population might be like.

##### Some terms
[1] Homocedasticity: A situation where the random disturbance between independent and dependent variables is the same across all values of the independent variables.

[2] Stationary Process: A stochastic process whose unconditional joint probability distributions between lags are the same.
