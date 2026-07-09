Numerical implementation of Prufer-based shooting for the Discrete Sturm-Liouville Eigenvalue Problem. 

Code and Supplemental materials for the paper: 

**About the paper**

"Theory and Computation of Discrete Sturm--Liouville Problems on Non-Uniform Grids via the Prufer Transformation"

Authors: Bailyn Hall(bhall76@ut.utm.edu, University of Tennessee at Martin), Kimsear Lor(klor@ut.utm.edu, University of Tennessee at Martin)

Project Advisors: Shalmali Bandyopadhyay (sbandyo5@utm.edu, University of Tennessee at Martin), Jacob Blazejewski (blazejewskijj@appstate.edu, Appalachian State University)

Program: National Research Experience Undergraduate Program, University of Tennessee at Martin, 2026

This repository contains Python implementations comparing Prufer-based shooting and regular shooting on a special case of the Sturm-Liouville eigenvalue problem: y″ = λy on [0, π].
Both methods are tested figures are generated for low eigenvalues(1-10), then higher eigenvalues(up to 60) on three grids: uniform grid, clustered grid, and graded grid.

**Requirements**

- Python 3.10+
- numpy
- scipy
- matplotlib

**Reproducing the results**

The Python files are divided by types of grids. Each Python file will include a block that tests regular shooting and another that tests Prufer-based shooting. Each block will generate their own data, format their own tables, and plot their own graphs. 

The **easiest** way to reproduce these results is to use Google Colab. By having a Google account, you have access to Google Colab. Copy a block of code that you desire to simulate and run it in the notebook. 

To test out other numbers of eigenvalues, look for num or numeig(current: either 10 or 60); that is where you input a value to get that desired number of eigenvalues. 

To change number of points on the graph, look for the variable n(default n = 100).

To change grid behaviour, i.e. number of points on the grids and spacing, modify the following variables for the following grids: 

- A variable for Uniform grid is number of points on the grid(n).
- Variables for Clustered grids are Fraction of Domain(frac_dom), Fraction of Points(frac_pts), and number of points on the grid(n). 
- Variables for graded grids are Ratio(ratio) and number of points on the grid(n).

Any other questions, contact Kimsear Lor(klor@ut.utm.edu). 

