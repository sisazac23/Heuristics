# Heuristics

This repository contains the codes developed for the Heuristics course at EAFIT University.

The purpose of the algorithms was to solve the travelling salesman problem for a specific data set proposed by the teacher.

TSP ask the following question: What is the shortest route that visits each city precisely once and returns to the starting city, given a list of cities and the distances between each pair of cities? It is an NP-hard combinatorial optimization problem. For the solution the course use the TSP integer linear programming formulation. Especifically, the Dantzig-Fulkerson-Johnson formulation.
Label the cities with the numbers $1, \ldots, n$ and define:

```math
x_{i j}= \begin{cases}1 & \text { the path goes from city } i \text { to city } j \\ 0 & \text { otherwise }\end{cases}
```

Take $c_{i j}>0$ to be the distance from city $i$ to city $j$. Then TSP can be written as the following integer linear programming problem:


$$\min \sum_{i=1}^n \sum_{j \neq i, j=1}^n  c_{i j} x_{i j}:$$

$$\sum_{i=1, i \neq j}^n x_{i j}=1   ~  ~ ~  ~  ~ ~  j=1, \ldots, n ;$$

$$\sum_{j=1, j \neq i}^n x_{i j}=1  ~  ~ ~  ~  ~ ~ i=1, \ldots, n$$

```math
\sum_{i \in Q} \sum_{j \neq i, j \in Q} x_{i j} \leq|Q|-1  ~  ~ ~  ~  ~ ~ \forall Q \subsetneq\{1, \ldots, n\},|Q| \geq 2
```
