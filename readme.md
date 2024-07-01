# TIF: Time-series-based Image Fusion

Welcome to the **TIF** repository! This package includes the source code and examples of the **Time-series-based Image Fusion (TIF)** algorithm. The TIF algorithm was developed to produce 10 m Harmonized Landsat and Sentinel-2 (HLS) data by fusing 30 m Landsat 8 and 10 m Sentinel-2 time series.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Overview

The TIF algorithm is designed to enhance the spatial resolution of Landsat 8 imagery from 30 m to 10 m by leveraging the higher resolution Sentinel-2 data. This results in 10 m harmonized data that is temporally consistent with both Landsat and Sentinel-2 observations.

## Features

- **Spatial Resolution Enhancement:** Improve Landsat 8 imagery resolution from 30 m to 10 m.
- **Temporal Consistency:** Harmonize data from different satellite sources to ensure temporal consistency.
- **Robust Performance:** Demonstrated robustness to temporal changes and varying land cover types.
- **Comprehensive Metrics:** Evaluation using RMSE, AAD, ERGAS, SSIM, and UIQI metrics to ensure high-quality results.

## Installation

To install the TIF package, you can clone the repository and ensure you have the required MATLAB toolboxes the [Mapping Toolbox](https://www.mathworks.com/products/mapping.html) and the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html).


```bash
git clone https://github.com/yourusername/TIF.git
```

## Usage
To use the TIF algorithm, follow these steps:

1. Prepare your Landsat 8 and Sentinel-2 time series data.
2. Use the TIF algorithm to fuse the data and produce harmonized 10 m HLS data.

Here's an example script to get you started:

## Examples
We have included several examples in the examples directory to demonstrate the usage of the TIF algorithm. These examples cover different scenarios and use cases, helping you understand how to apply the algorithm to your own data.

### Example 1: Basic Usage 
This example demonstrates the basic usage of the TIF algorithm with a single pixel.


### Example 2: Advanced Usage
This example shows advanced usage with additional options and parameters.


## License
## Acknowledgement




