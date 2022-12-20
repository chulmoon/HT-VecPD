# Hypothesis Testing for Vectorized Persistence Diagram
Topological data analysis involves the statistical characterization of the shape of data. Persistent homology is a primary tool of topological data analysis, which can be used to analyze topological features and perform statistical inference. We present a two-stage hypothesis test for vectorized persistence diagrams. The first stage filters vector elements in the vectorized persistence diagrams to enhance the power of the test. The second stage consists of multiple hypothesis tests, with false positives controlled by false discovery rates. We demonstrate the flexibility of our method by applying it to a variety of simulated and real-world data types. Our results show that the proposed hypothesis test enables accurate and informative inferences on the shape of data compared to the existing hypothesis testing methods for persistent homology.

# Simulation 
## “./simulation_point_cloud” folder
* `data_generation.R`: Generate point cloud data and compute persistence diagrams
* `data_persistence_image.R`: Generate persistence images
* `functions.R`: Functions for hypothesis tests
* `test_kernel.ipynb`: Conduct hypothesis tests using persistence weighted Gaussian kernel (PWGK)
* `test_persistence_diagram.R`: Conduct hypothesis tests using persistence diagrams
* `test_persistence_landscape.R `: Conduct hypothesis tests using persistence landscapes
* `test_two_stage.R `: Conduct hypothesis tests using two-stage hypothesis testing

## “./simulation_beetle” folder
* `data_persistence_diagram.R`: Generate beetle population data and compute persistence diagrams
* `functions.R`: Functions for hypothesis tests
* `test_dtw.R`: Conduct hypothesis tests using dynamic time warping (DTW)
* `test_kernel.ipynb`: Conduct hypothesis tests using persistence weighted Gaussian kernel (PWGK)
* `test_persistence_diagram.R`: Conduct hypothesis tests using persistence diagrams
* `test_persistence_landscape.R `: Conduct hypothesis tests using persistence landscapes
* `test_two_stage.R `: Conduct hypothesis tests using two-stage hypothesis testing

## “./simulation_rock” folder
* `rock_simulation_datagen.ipynb`: Generate binary rock image data and compute persistence diagrams
* `data_persistence_diagram.R`: Process persistence diagram
* `functions.R`: Functions for hypothesis tests
* `test_kernel.ipynb`: Conduct hypothesis tests using persistence weighted Gaussian kernel (PWGK)
* `test_persistence_diagram.R`: Conduct hypothesis tests using persistence diagrams
* `test_persistence_landscape.R `: Conduct hypothesis tests using persistence landscapes
* `test_two_stage.R `: Conduct hypothesis tests using two-stage hypothesis testing

# Application
## “./application_rock” folder
* `data_persistence_diagram.ipynb`: Compute persistence diagrams
* `data_persistence_image.R`: Convert persistence diagrams to persistence images
* `functions.R`: Functions for hypothesis tests
* `permulabel.R`: Permutation for hypothesis tests using persistence diagrams
* `test_kernel.ipynb`: Conduct hypothesis tests using persistence weighted Gaussian kernel (PWGK)
* `test_persistence_diagram.R`: Conduct hypothesis tests using persistence diagrams
* `test_persistence_landscape.R `: Conduct hypothesis tests using persistence landscapes
* `test_two_stage.R `: Conduct hypothesis tests using two-stage hypothesis testing

## “./application_sound” folder
* `functions.R`: Functions for hypothesis tests
* `test_two_stage.R `: Conduct hypothesis tests using two-stage hypothesis testing
