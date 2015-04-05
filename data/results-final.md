# Test results

All results are evaluated according to the mean miss rate on the test set.

#### De Falco et al. results (expected)

| Database  | Ψ1  | Ψ2 | Ψ3 |
| :------------ |:------:| :-----:|  :-----:|
| iris          | 3.68% | **2.63%** | 5.26% |
| wine          | 6.22% | **2.22%** | 2.88% |
| glass         | 40.18% | 39.05% | **38.67%** |
| thyroid       | 6.66% | 5.55% | **3.88%** |

#### Measured results

| Database  | Ψ1 | Ψ2 | Ψ3 |
| :------------ |:------:| :-----:|  :-----:|
| iris          | 6.02% | 7.94% | **4.74%** |
| wine          | 12.88% | **2.55%** | 2.89% |
| glass         | TBT | TBT | TBT |
| thyroid       | TBT | TBT | TBT |

> *TBT: To be tested.*

###### Test details

**Wine**

| fitnessfcn  | stop_cause | std | gens |
| :------------ |:------:| :-----:|  :-----:|
| Ψ1    | 1 | ±5.88 | 239 |
| Ψ2    | 2 | ±1.94 | 355 |
| Ψ3    | 2 | ±2.12 | 387 |

> Stop causes:<br>
> 1: The value of the fitness function did not improve in the last *n* generations.<br>
> 2: Average cumulative change in value of the fitness function over *n* generations less than *x*.<br>
> 3: Reached limit of *n* iterations.<br>