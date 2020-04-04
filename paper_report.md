
# Analysis and Report on the Paper "Energy Management Optimization in a Battery/Supercapacitor Hybrid Energy Storage System"

M. Choi, S. Kim and S. Seo, "Energy Management Optimization in a Battery/Supercapacitor Hybrid Energy Storage System," in IEEE Transactions on Smart Grid, vol. 3, no. 1, pp. 463-472, March 2012.

DOI: 10.1109/TSG.2011.2164816

## Issues / Mistakes in the Paper


It seems there might be some issues in the paper. At this point it is not sure how they will influence the results/calculation. The issues are as follows:

- The value for Delta is given nowhere in the paper besides the info "We assume that the voltage and current are measured using discrete signals with a sufficiently small sampling period". The stepping for $t \in T$ is $1$ therefore I set $\Delta$ to `Delta = 200e-6$;` (200 microseconds), since that seemed a reasonable number for a sampling interval for currents.

- The "equivalent cost function" mentioned in the paper $f(\cdot)$ for finding $\epsilon$ is not defined. However they define a value for epsilon they used. This value is $\epsilon=0.7$ which I used.

- For the MIAD algorithm (Algorithm 2 in the paper), the initial values for $temp_14 and $temp_2$ are not set / not defined. I set them to be $\sigma_1$ and $\sigma_2$ initially, since they will be overwritten after the first iteration within the MIAD algorithm anyways.

- Within the proof of 2nd constraint in the appendix in the fifth line, the equation seems to be missing a minus sign between the two terms $I_{s_k}^{out}$ and $I_{s_k}^{in}$. I think the author just missed the minus while typing the equation out. Because otherwise the equation would not be equal to the equation in the fourth line

- Within the proof of 2nd constraint in the appendix the author eliminates the $|\ldots|$ function between line 3 and and 4. The author switches the sign between the two terms from negativ to positive. I assume to presever the nature of the problem (resistive energy losses are always > 0). However this is only valid as long as $I_{s_k}^{out} \leq I_{s_k}^{in}$, which is not incorporated anywhere in the paper. I used the equations from the paper whithout any additional inequality condition.

- The authors give no parameters which they used for the solver they used. Therefore I am in the progress of exploring parameters now.

## Logic / Algorithm Implmeneted in the Paper

Now what I think about the logic in the paper:
- The authors minimize a pretty complex objective/penalty function, they do this because they claim it reflects how the battery would deteriorate with RMS current and fluctuations. From my experience it does not matter, which objective function to minimize, as long as the objective function represents the goals, therefore  I think using a simpler objective function would yield the same results but make computation easier.

- The authors introduced $L_k$ in **P2** (and therefore also **P3**) in order to map the **P2** problem to a LP optimization problem. As I see it: when combining **P1** and **P2** to **P3**, **P1** is not LP but a convex problem, therefore we gain nothing by mapping **P2** to a LP problem, and it is not necessary to introduce the additional L_k terms, which result in the additional inequality constraint $A\cdot x <= b$ as well as $K\cdot T$ more variables in $x$. Instead it should be possible to directly minimize the $sum_{k\in K} || I_sk_out - I_sk_in ||_1 $ term.

## Input Values
The authors do not give any values for M1-M6, only plots of those. 
I reverse engineered M1-M3 and came up with some reasonable arrays for M4-M6, which look similar to the ones used in the paper.

Also I used symbolic math in python to check/proof that the Algorithms I implemented reflect those given in the paper. This is given in the jupyter notebook.