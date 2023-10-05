#!/usr/bin/env python3

import svd as mysvd
import numpy as np

def get_masses():
    """
    Get the number of masses from the user.
    
    Returns:
        int: Number of masses.
    """
    while True:
        try:
            num_masses = int(input("Enter the number of masses: "))
            if num_masses > 0:
                break
            else:
                print("Number of masses must be a positive integer. Please try again.")
        except ValueError:
            print(
                "Invalid input. Please enter a positive integer for the number of masses.")
    return num_masses

def get_springs():
    """
    Get the number of springs from the user.
    
    Returns:
        int: Number of springs.
    """
    while True:
        try:
            num_springs = int(input("Enter the number of springs: "))
            if num_springs > 0:
                break
            else:
                print("Number of springs must be a positive integer. Please try again.")
        except ValueError:
            print(
                "Invalid input. Please enter a positive integer for the number of springs.")
    return num_springs
    
def get_input():
    """
    Collect  inputs from the user for the spring-mass system.
    
    Returns:
        num_springs (int): the number of springs
        num_masses (int): the number of masses
        spring_constants (list): list of the spring constants
        masses (list): list of the masses (kg)
        boundary condition (int): the boundary condition where '1' is fixed/fixed and '2' is fixed/free
    """

    # Get the number of springs
    num_springs = get_springs()
    num_masses = get_masses()

    # Get the number of masses


    # Ensure that the number of springs and masses are compatible
    while num_springs != num_masses and num_springs != num_masses + 1:
        print("The number of springs and masses are not compatible. Ensure that:")
        print("- The number of springs is equal to the number of masses, or")
        print("- The number of springs is equal to the number of masses plus one.")
        num_springs = get_springs()
        num_masses = get_masses()
        

    # Get the spring constants
    spring_constants = []
    for i in range(num_springs):
        while True:
            try:
                k = float(input(f"Enter the spring constant for spring {i+1}: "))
                if k > 0:
                    spring_constants.append(k)
                    break
                else:
                    print("Spring constant must be a positive number. Please try again.")
            except ValueError:
                print("Invalid input. Please enter a positive number for the spring constant.")

    # Get the masses
    masses = []
    for i in range(num_masses):
        while True:
            try:
                m = float(input(f"Enter the mass for mass {i+1}: "))
                if m > 0:
                    masses.append(m)
                    break
                else:
                    print("Mass must be a positive number. Please try again.")
            except ValueError:
                print("Invalid input. Please enter a positive number for the mass.")

    # Get the type of boundary conditions
    while True:
        bc = input("Enter the type of boundary condition (1 or 2 fixed ends):\nNote: 1 means 'fixed/fixed' system, 2 means 'fixed/free' system\n")
        if int(bc) in [1, 2]:
            if ((int(bc) == 1) and (int(num_masses) != int(num_springs) - 1)) or ((int(bc) == 2) and (int(num_masses) != int(num_springs))):
                print(
                    "Invalid input. Springs and Masses do not match BC. Please enter 1 or 2.")
            else:
                break
        else:
            print("Invalid input. Please enter 1 or 2.")

    return int(num_springs), int(num_masses), spring_constants, masses, bc

def get_A(num_springs, num_masses):
    """
    Generate a matrix A based on the number of springs and masses.
    
    Parameters:
        num_springs (int): Number of springs.
        num_masses (int): Number of masses.
        
    Returns:
        np.ndarray: the A matrix (the difference matrix).
    """
    # Initialize A with zeros
    A = np.zeros((num_springs, num_masses))

    # Insert 1 to the main diagonal
    np.fill_diagonal(A, 1)

    # Insert -1 to the lower diagonal if number of springs is more than 1
    if num_springs > 1:
        np.fill_diagonal(A[1:], -1)

    return A

def get_C(spring_constants):
    """
    Create a diagonal matrix with spring constants.
    
    Parameters:
        spring_constants (list): List of spring constants.
    
    Returns:
        np.ndarray: Diagonal matrix with spring constants.
    """
    return np.diag(spring_constants)

def get_stiff(A, C):
    """
    Calculate the stiffness matrix (K) using matrices A and C.
    
    Parameters:
        A (np.ndarray): Matrix derived from number of springs and masses.
        C (np.ndarray): Diagonal matrix of spring constants.
        
    Returns:
        K (np.ndarray): Stiffness matrix.
    """
    return (A.T @ C @ A)

def get_F(masses):
    """
    Calculate force due to gravity on each mass.
    
    Parameters:
        masses (list): List of masses.
    
    Returns:
        F (np.ndarray): Force vector.
    """
    g = 9.81  # gravity in m/s^2
    F = np.array(masses) * g
    return F

def get_displacements(K, f):
    """
    Calculate displacements for given stiffness matrix and force vector.
    
    Parameters:
        K (np.ndarray): Stiffness matrix.
        f (np.ndarray): Force vector.
    
    Returns:
        U (np.ndarray): Displacement vector.
    """
    K_inv = mysvd.get_Ainv(K)
    U = K_inv @ f
    return U

def get_elongations(A, u):
    """
    Calculate elongations of springs using matrix A and displacement vector u.
    
    Parameters:
        A (np.ndarray): Matrix derived from the number of springs and masses.
        u (np.ndarray): Displacement vector.
    
    Returns:
        e (np.ndarray): Elongation vector.
    """
    return A @ u

def get_internalStresses(C, e):
    """
    Calculate internal stresses in the springs using matrix C and elongation vector e.
    
    Parameters:
        C (np.ndarray): Diagonal matrix of spring constants.
        e (np.ndarray): Elongation vector.
    
    Returns:
        w (np.ndarray): Internal stress vector.
    """
    return C @ e

def solve_system(num_springs, num_masses, spring_constants, masses, bc):
    """
    Solve the spring-mass system given the number of springs, masses, spring constants, 
    masses, and boundary conditions. This function computes several important matrices 
    and vectors used in the analysis of the spring-mass system, such as the difference 
    matrix, spring constant matrix, stiffness matrix, force vector, displacement vector,
    elongation vector, and internal stress vector.
    
    Parameters:
        num_springs (int): The number of springs in the system.
        num_masses (int): The number of masses in the system.
        spring_constants (list of float): The spring constant for each spring in the system.
        masses (list of float): The mass value for each mass in the system.
        bc (int): An integer representing the type of boundary condition (1 for fixed/fixed, 2 for fixed/free).
    
    Returns:
        tuple: A tuple containing the following elements in order:
            - A (np.ndarray): The difference matrix.
            - C (np.ndarray): The diagonal spring constant matrix.
            - K (np.ndarray): The stiffness matrix.
            - F (np.ndarray): The force vector.
            - u (np.ndarray): The displacement vector.
            - e (np.ndarray): The elongation vector.
            - w (np.ndarray): The internal stress vector.
            
    Note:
        Ensure that the input lists `spring_constants` and `masses` have lengths equal 
        to `num_springs` and `num_masses` respectively, and that the elements are non-negative.
        Boundary condition, `bc`, should be either 1 or 2, otherwise, it may not behave as expected.
    """
    A = get_A(num_springs, num_masses)
    C = get_C(spring_constants)
    K = get_stiff(A, C)
    F = get_F(masses)
    u = get_displacements(K, F)
    e = get_elongations(A, u)
    w = get_internalStresses(C, e)
    return A, C, K, F, u, e, w


def main():

    # could just do:
    # A, C, K, F, disp, elong, inStress = solve_system(get_input())

    # Creating system
    num_springs, num_masses, spring_constants, masses, bc = get_input()

    # Printing out system
    print('-'*40)
    print(f'Number of Springs: {num_springs}')
    print(f'Number of Masses: {num_masses}')
    print(f'Spring Constants: {spring_constants}')
    print(f'Masses: {masses}')
    print(f'Boundary Condition: {bc}')


    # Calculating System
    A, C, K, F, disp, elong, inStress = solve_system(num_springs, num_masses, spring_constants, masses, bc)

    # # # Printing out Solutions
    # print(mysvd.pretty_matrix(A, "A"))
    # print(mysvd.pretty_matrix(C, "C"))
    # print(mysvd.pretty_matrix(K, "K"))
    # print(mysvd.pretty_matrix(F, "F"))
    mysvd.pretty_matrix(disp, "Equilibrium Displacement")
    mysvd.pretty_matrix(inStress, "Internal Stress")
    mysvd.pretty_matrix(elong, "Elongation")

    cond_num = mysvd.get_condition(K) # This is the function I made
    print('-'*40, f'\nCondition Number: ', ('%.5f' % cond_num), '\n', '-'*40)

if __name__ == "__main__":
    main()
