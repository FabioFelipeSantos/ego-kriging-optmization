<style>
    * {
        line-height: 2.2
    }
    code {
        padding: 3px 8px
    }
</style>

# Tutorial for Using the Optimization Software with EGO and Kriging

This tutorial is a step-by-step guide on how to use the **Efficient Global Optimization (EGO)** software that utilizes surrogate models with Kriging.

A large part of this optimizer's development was based on the following texts:

-   FORRESTER, A.; SOBESTER, A. S.; KEANE, A. _Engineering design via surrogate modelling: a practical guide_. Chichester, West Sussex, United Kingdom: John Wiley & Sons, 2008.
-   JONES, D. R. A taxonomy of global optimization methods based on response surfaces. _Journal of Global Optimization_, v. 21, n. 4, p. 345-383, 2001.
-   JONES, D. R.; SCHONLAU, M.; WELCH, W. J. Efficient global optimization of expensive black-box functions. _Journal of Global Optimization_, v. 13, n. 4, p. 455-492, 1998.

For any errors or suggestions, please contact the creators:

1.  **Fábio Felipe dos Santos**: Graduate Program in Mechanical Engineering, POSMEC UFSC, fabiotrabmat@gmail.com;
2.  **Dr. Rafael Holdorf Lopez**: Department of Civil Engineering and Graduate Program in Civil Engineering, PPGEC UFSC, rafaelholdorf@gmail.com.

The tutorial is divided into two parts. The first part covers the adjustment of the optimizer and objective function parameters. In the second part is provided some advices about the optimization proccess and results.

## PART 1

There is in this repository root a folder named `EGO-Kriging`. This is main folder were is the software and it contains the following:

-   The file `EGO_Kriging.m` responsible for defining the user's optimization problem and parameters.
-   The `script_opt` subfolder contains all the optimizer routines in `.m` files.
-   The `fobj` subfolder should store the user's objective functions in `.m` files.

### 1. Parameter Adjustment

First, it is recommended to edit the optimization parameters inside the `EGO_Kriging.m` file.

Those parameters must be passed as a `struct` to the optimization routine, so it is recommended to use the format `StructName.paramName = paramValue`.

The parameter names must be exactly the same as those listed in Tables 1.1 and 1.2.

#### Table 1.1: Mandatory Parameters

| Name                | Description                                                                                    |
| :------------------ | :--------------------------------------------------------------------------------------------- |
| `func`              | Name of the objective function (must be a string).                                             |
| `maxFE`             | Maximum number of objective function evaluations.                                              |
| `k`                 | Number of dimensions of the problem.                                                           |
| `LB`                | A $1 \times k$ vector with the lower bounds of the variables.                                  |
| `UB`                | A $1 \times k$ vector with the upper bounds of the variables.                                  |
| `calculoVetorizado` | Defines if the objective function can be calculated in vector form. More details in Section 2. |

#### Table 1.2: Optional Parameters

| Name         | Description                                                                                                                                                            |
| :----------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `n`          | Number of points in the initial sample. If not provided, the following rule will be used by default: $n=6.5k$ if $k \le 5$ or $n=4.5k$ if $k \gt 5$, rounded up.       |
| `S`          | Initial sample with `n` points belonging to $[0,1]^k$. Default: A set determined by a Latin Hypercube.                                                                 |
| `tolMinEI`   | Minimum allowed tolerance for the Expected Improvement. Default value: $10^{-40}$.                                                                                     |
| `nPop`$^{*}$ | Number of individuals for the GA heuristic optimizer. Default value: $20 \cdot k$.                                                                                     |
| `graf`       | Creates a graphical output for 1D and 2D functions. Only two possible values: `0` (default) - No graphical output . `1` - Displays the graphical output on the screen. |

\* **Note:** MATLAB's own GA optimizer is used in two different and extremely essential processes: determining the surrogate model parameters and determining the Infill Point (IP). A too lower value will generate a bad surrogate and a incorrect infill point in EGO. A too upper value will considerably increase the computational time.

The optimizer is called using the `otimEGO.m` function (the last line of code in `EGO_Kriging.m`).

After finished, the software should return a variable called `Resultado` being a struct with the following variables:

#### Table 1.3: Description of Output Variables in `Resultado`

| Name                | Description                                             |
| :------------------ | :------------------------------------------------------ |
| `MelhorMinimizador` | Best minimizer found.                                   |
| `MelhorValorOtimo`  | Best minimum value found.                               |
| `PosMelhorMin`      | Iteration in which the best minimum value was obtained. |
| `numFE`             | Number of objective function evaluations.               |
| `numIPAdd`          | Number of Infill Points that were added.                |
| `numPontosIniciais` | Number of points in the initial sample.                 |
| `fobj`              | Contains the objective function values at all points.   |
| `Histfmin`          | History of the best minimums.                           |
| `Histdmin`          | History of the minimizers.                              |
| `CriterioParada`    | Reason for the optimizer's termination.                 |

The `fitKrg` variable is a `struct` that contains all the necessary information about the final surrogate model obtained.

### 2. Objective Function

The objective function must be specified in a `.m` type file.

This function must receive a single argument and return only a single argument.

```matlab
function f = FunctionName(d)
    % Implement the code
end
```

where `d` is the input argument and `f` is the output argument.

If the `calculoVetorizado` parameter has the value `0`, then the function needs to handle a single vector of design variables, returning only a single value for the objective function. Otherwise, with `calculoVetorizado` setted to value `1`, then your function must be implemented with the vectorized evaluation of a $n \times k$ matrix $d$ of $n$ design inputs with $k$ design variables. The returned value must be a column vector with $n$ rows.

As an example, suppose we want to optimize the function $f(d) = d_{1}^{2} + d_{2}^{2} + d_{3}^{2}$, where $d = \{d_1, d_2, d_3\}$ is a vector of 3 design variables.

For the non-vectorized case (`calculoVetorizado = 0`), we suggest:

```matlab
function f = FunctionName(d)
    f = d(1)^2 + d(2)^2 + d(3)^2;
end
```

For the vectorized case (`calculoVetorizado = 1`), we suggest:

```matlab
function f = FunctionName(d)
    f = d(:,1).^2 + d(:,2).^2 + d(:,3).^2;
end
```

---

## PART 2

This part contains some tips regarding the parameters and the results that are expected from the optimizer.

-   The optimizer was designed to work well as the number of dimensions in the design variable is under 20. Otherwise, the computational time can be extremely prohibitive

-   The number of initial points `n` does not have a fixed value, but it should be appropriate for the complexity of the problem.

-   The initial sample `S`, if not declared, is obtained from a Latin Hypercube within the $[0,1]^k$ hyperspace.

-   The number of added IPs is crucial for obtaining the best minimum.

-   After the optimization is finished, the user can make new predictions using the Kriging model based on the information from the last created surrogate model. This can be done by using the following call:

    <div style="widht: 100%; text-align: center; margin-bottom: 6px"><code style="padding: 6px 12px">[predicao, erroPred] = predKrg(d, fitKrg)</code></div>

    where `predicao` is the Kriging prediction value and `erroPred` is the variance of the prediction error, both calculated at point `d`, where $\bold{d} \in [0,1]^k$.

-   To calculate the Expected Improvement, one can use the call:

    <div style="widht: 100%; text-align: center; margin-bottom: 6px"><code style="padding: 6px 12px">EI = -expImp(d, fitKrg)</code></div>

    where the negative sign is necessary because the Expected Improvement value is multiplied by -1 for its maximization by the GA.

-   Since the entire optimization process is based on random numbers, especially in the initial population and the GA, it is very likely that the minimum value obtained will be different each time the algorithm is run.
