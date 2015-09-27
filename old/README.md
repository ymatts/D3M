# D3M
R source code of D3M.</br>
Authour: Yusuke Matsui</br>
Email: ymatsui@med.nagoya-u.ac.jp</br>

#Package information:
Dependency: MASS</br>

This directory is old version of D3M.
This package consists of four files.</br>
#1. sample.txt: </br>
Simulation dataset</br>
Data is 8 * 300 matrix: The rows correspond to 8 cases. The columns correspond to 160, 140 samples of case and control group respectively.</br>
#2. d3m.R:</br>
Main code of d3m. When using it, please include the source code by 'source("d3m.R")'.</br>
#3. plot_d3m.R: </br>
The code for visualizing distribution of beta values with Beanplot and Q-Q plot. When using it, please include the source code by 'source("plot_d3m.R")'.</br>
#4. example.R: </br>
A sample code demonsrating the functions above.</br>

