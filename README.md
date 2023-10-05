# SVD Spring Mass System
This repository contains python 2 python scripts to preform singular value decomposition (SVD) on a matrix and to solve relevent values of a spring-mass system. 

## Scripts
### `svd.py`
This script prefors SVD on a given matrix. Finding, A = U &Sigma; V<sup>T</sup> 

The SVD method decomposes a matrix into three other matrices: U, Sigma, and V. U represents the left singular vector matrix, V<sup>T</sup> is the right singular vector matrix, and Sigma is a diagonal matrix containing the singular values, or the eigenvalues of AA<sup>T</sup> and A<sup>T</sup>A, which are equal. This script also compares the found SVD values to the `numpy.linalg.svd()` blackbox function.


The 3 most notable functions from this script are:

`svd(A)`: Perform Singular Value Decomposition on matrix A; returning U, Sigma, and V.

`get_condition(A)`: Compute the condition number of matrix A using its singular values to find the spectral norm of the matrix. 

`get_Ainv(A)`: Calculate the inverse of matrix A, if it exists. (If a singular value = 0, the inverse does not exist)

This code is written to take user input when executed directly.

**Example (for a non-singular matrix):**

```
Jacksons-MBP:project_1 jacksonthetford$ python3 svd.py
Enter the number of rows: 2
Enter the number of columns: 3
Enter element at position 1, 1: 3
Enter element at position 1, 2: 2
Enter element at position 1, 3: 2
Enter element at position 2, 1: 2
Enter element at position 2, 2: 3
Enter element at position 2, 3: -2
---------------------------------------- 
Input Matrix:

[[ 3.  2.  2.]
 [ 2.  3. -2.]]

---------------------------------------- 
U Matrix:

[[-0.70711  0.70711]
 [-0.70711 -0.70711]]

---------------------------------------- 
Sigma Matrix:

[[5. 0. 0.]
 [0. 3. 0.]]

---------------------------------------- 
V.T Matrix:

[[-0.70711 -0.70711 -0.     ]
 [ 0.2357  -0.2357   0.94281]
 [-0.66667  0.66667  0.33333]]

---------------------------------------- 
U @ Sigma @ V.T (should be A) Matrix:

[[ 3.  2.  2.]
 [ 2.  3. -2.]]

---------------------------------------- 
Condition Number:  1.66667
---------------------------------------- 
A Inverse Matrix:

[[0. 0.]
 [0. 0.]
 [0. 0.]]

---------------------------------------- 
Blackbox Input Matrix:

[[ 3.  2.  2.]
 [ 2.  3. -2.]]

---------------------------------------- 
Blackbox U Matrix:

[[-0.70711 -0.70711]
 [-0.70711  0.70711]]

---------------------------------------- 
Blackbox Sigma Vector:

[5. 3.]

---------------------------------------- 
Blackbox V Matrix:

[[-0.70711 -0.70711 -0.     ]
 [-0.2357   0.2357  -0.94281]
 [-0.66667  0.66667  0.33333]]
```

**Example comparing my code to the blackbox for an invertible matrix:**

Using the matrix, `A = [[1, 0], [1, 0]]`

```
Jacksons-MBP:project_1 jacksonthetford$ python3 svd.py
Traceback (most recent call last):
  File ".../coe352/project_1/svd.py", line 138, in get_Ainv
    raise Exception(
Exception: Error: Input matrix is singular and cannot be inverted
Jacksons-MBP:project_1 jacksonthetford$ python3 svd.py
---------------------------------------- 
Blackbox Input Matrix:

[[1 0]
 [1 0]]

---------------------------------------- 
Blackbox U Matrix:

[[-0.70711 -0.70711]
 [-0.70711  0.70711]]

---------------------------------------- 
Blackbox Sigma Vector:

[1.41421 0.     ]

---------------------------------------- 
Blackbox V Matrix:

[[-1. -0.]
 [ 0.  1.]]

Jacksons-MBP:project_1 jacksonthetford$ 
```
My svd function will catch the singular nature of the matrix and raise an exception, however the `numpy.linalg.svd` function will not, and will still provide an answer that is incorrect. 


### `spring_mass.py`
This script accepts user input for a spring mass system and then solves the system by calculating the equilibrium displacements, the internal stresses and the elongations of the spring/mass system. To solve this system, the svd function from my svd.py is used to invert the stiffness matrix to solve **f = Ku**. Where **f** is the force vector, **K** is the stiffness matrix, and **u** is the equilibrium displacements. 


The code allows for user input for the number of springs/masses, the spring constants for each spring, the masses, and which boundary condition to apply (either fixed at both ends or just one). 

This script solves the system by the following method, using the `solve_system` function which takes the user input: `num_springs`, `num_masses`, `spring_constants`, `masses`, `bc`.

First, the A, or difference, matrix is found as it represents the relationship between the displacments and elongations of the spring as **e = A u**. Then, the C matrix is found; it is a diagonal matrix of the user-inputted spring constants. The K matrix is then found through the identity, K = A<sup>T</sup> C A. The f vector, the gravitational force on the masses, is the m<sub>i</sub>*g, given the user-input for the masses and the `svd.get_Ainv` function is utilized to find the equilibrium displacement, u, such that as **f = K u**, u can be found by **K<sup>-1</sup> f = u**. From this, u is known so the elongations can be found (from **e = A u**), then the internal stresses can be found by **w = C e**.

An example of running this script using a 4 spring, 3 mass, fixed/fixed system:

```
Jacksons-MBP:project_1 jacksonthetford$ python3 spring_mass.py 
Enter the number of springs: 4
Enter the number of masses: 3
Enter the spring constant for spring 1: 1
Enter the spring constant for spring 2: 2
Enter the spring constant for spring 3: 3
Enter the spring constant for spring 4: 4
Enter the mass for mass 1: 10
Enter the mass for mass 2: 10
Enter the mass for mass 3: 20
Enter the type of boundary condition (1 or 2 fixed ends):
Note: 1 means 'fixed/fixed' system, 2 means 'fixed/free' system
1
----------------------------------------
Number of Springs: 4
Number of Masses: 3
Spring Constants: [1.0, 2.0, 3.0, 4.0]
Masses: [10.0, 10.0, 20.0]
Boundary Condition: 1
---------------------------------------- 
A Matrix:

[[ 1.  0.  0.]
 [-1.  1.  0.]
 [ 0. -1.  1.]
 [ 0.  0. -1.]]

None
---------------------------------------- 
C Matrix:

[[1. 0. 0. 0.]
 [0. 2. 0. 0.]
 [0. 0. 3. 0.]
 [0. 0. 0. 4.]]

None
---------------------------------------- 
K Matrix:

[[ 3. -2.  0.]
 [-2.  5. -3.]
 [ 0. -3.  7.]]

None
---------------------------------------- 
F Vector:

[ 98.1  98.1 196.2]

None
---------------------------------------- 
Equilibrium Displacement Vector:

[0. 0. 0.]

---------------------------------------- 
Internal Stress Vector:

[0. 0. 0. 0.]

---------------------------------------- 
Elongation Vector:

[0. 0. 0. 0.]

---------------------------------------- 
Condition Number: 7.75840 
 ----------------------------------------
```

It should be notted that when the number of springs is 1 more than the number of masses, the system is always fixed at both ends. When the number of springs is equal to the number of masses, the system is always free at one end. No other configuration of springs and masses is physically possible. **In the case where both ends are free**, for out case of a vertical system, this would mean that the system at a whole would be in free-fall and that in a vacuum the internal forces would be zero and there would be no displacement, elongation, or internal stresses. In other words, the system is numerically indeterminate as there is simply not enough equations of equilibrium to determine the reactions at the supports. Also, this situation is impossible to achieve as there were always be some external force acting on the system (normal force, pressure, etc.) unless one went to space and let the spring mass float around.


#### More examples:
Note: The user input and the A, C, K matricies and the F vector output have been ommited.
##### Fixed/Free:













