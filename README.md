# Battery Current Flow Optimization using MATLAB
This software tries to implement the battery current flow optimization scheme shown in
https://ieeexplore.ieee.org/document/6046111/. 
However the results obtained did not match the results shown in the paper.

The main scripts are:
- `optim_P3.m`: for optimizing the problem P3 as given in the paper.
- `optim_P3_small.m`: small version of aforementioned script where all functions are packed in one file.
- `optim_P3_adj.m`: an alternative math implementation for solving the problem by me loosly based on the paper.

The math was checked 2 times with symbolic math toolboxes in python (file: `math_check_paper_2.ipnb`) and MATLAB (files: `SymCheck.m` and the outputs `symbolic_math_check_output.txt` as well as `html/SymCheck.html`) and it was found that the math implemented in this software matches the equations given in the paper. 

Whether the problem lies within a misinterpretation of the informaton given in the paper by me, or due to different initial conditions and parameters can not be said at this point. I am making this repository open source in the hope, that it will help people in future. 

A report of my findings when analyzing the paper is given in `paper_report.md` I do not guarantee for correctness though.

When checking the math in the paper before I got a suspicion about a slight error in the math. It is regarding the Appendix where they generate the constraint equation. There they replaced the `abs( I_out - I_in )` function with `( I_out + I_in )`. At this point they already limited `I_out >= 0` and `I_in >= 0`. However I have the suspicion that this step is not valid. But I doubt that this is the reason for the discrepancy in the results given. 

Another thing I noticed is:
`vk0` (voltage now) is set as a free variable for the optimization to make the math work. This however is not possible since one can not influence the voltage now anymore.
