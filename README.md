# D3M
R source code of D3M.
Authour: Yusuke Matsui
Email: ymatsui@med.nagoya-u.ac.jp

#Package information:
Dependency: MASS
This package consists of four files.
1. sample.txt: Simulation dataset
	Data is 8 * 300 matrix: The rows correspond to 8 cases in Table 1 of section 3. The columns correspond to 160, 140 samples of case and control group respectively.
2. d3m.R: Main code of d3m. When using it, please include the source code by 'source("d3m.R")'.
3. plot_d3m.R: The code for visualizing distribution of beta values with Beanplot and Q-Q plot. When using it, please include the source code by 'source("plot_d3m.R")'.
4. example.R: A sample code demonsrating the functions above.

