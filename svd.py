#!/usr/bin/env python3

import numpy as np
from math import sqrt

def input_matrix():
    '''
    Retrieve matrix details from user input.

    This function prompts the user to input the number of rows and columns 
    and each element of the matrix sequentially. It enforces input validity 
    checks to ensure all inputs are integers and floats for the dimensions 
    and elements, respectively.

    Returns:
        A (np.array): The user-defined matrix.
    '''
    global rows
    global cols
    while True:
        try:
            rows = int(input("Enter the number of rows: "))
            cols = int(input("Enter the number of columns: "))
            if rows <= 0 or cols <= 0:
                print("Please enter positive integer values for rows and columns.")
                continue
            break
        except ValueError:
            print("Please enter valid integer values.")
            continue

    matrix = []
    for i in range(rows):
        row = []
        for j in range(cols):
            while True:
                try:
                    element = float(
                        input(f"Enter element at position {i+1}, {j+1}: "))
                    break
                except ValueError:
                    print("Please enter a valid number.")
                    continue
            row.append(element)
        matrix.append(row)
    return np.array(matrix)

def svd(A: np.array) -> tuple:
    '''
    Perform Singular Value Decomposition on matrix A.

    This function decomposes matrix A into its singular value components
    U, Sigma, and V using eigenvalue decomposition and square root operation.

    Parameters:
        A : np.array
            The input matrix to decompose.

    Returns:
        tuple (U, Sigma, V)
            U : np.array
                Left singular vector matrix.
            Sigma : np.array
                Diagonal matrix containing singular values.
            V : np.array
                Right singular vector matrix.
    '''
    rows = A.shape[0]
    cols = A.shape[1]
    # Getting eigenvalues and eigenvectors
    eigVal_U, eigVec_U = np.linalg.eig(A@A.T)
    eigVal_V, eigVec_V = np.linalg.eig(A.T@A)

    # Getting Sigma
    sigma = np.zeros((rows, cols))
    for i in range(min(rows, cols)):
        sigma[i][i] = sqrt(eigVal_U[i])

    # Sort the eigen values 
    sort_index_U = eigVal_U.argsort()[::-1]
    eigVal_U = eigVal_U[sort_index_U]
    eigVal_U = eigVec_U[:, sort_index_U]

    sort_index_V = eigVal_V.argsort()[::-1]
    eigVal_V = eigVal_V[sort_index_V]
    eigVec_V = eigVec_V[:, sort_index_V]

    # Get U and V (already sorted), then correct sign of U
    U = eigVec_U
    V = eigVec_V
    same_sign = np.sign((A @ V)[0][0] * (U @ (eigVal_U @ np.eye(rows)))[0][0])
    U = U * same_sign.reshape(1, -1)

    return U, sigma, V

def get_condition(A: np.array):
    '''
    Compute the condition number of matrix A using its singular values.

    This function calculates the condition number of A by obtaining the 
    singular values through eigenvalue decomposition and returns their ratio 
    (max/min).

    Parameters:
        A (np.array): Input matrix for which to calculate the condition number.

    Returns:
        cond_num (float): Condition number of matrix A.
    '''
    eigVal_U, _ = np.linalg.eig(A@A.T)
    sigma = np.sqrt(eigVal_U)

    return (max(sigma)/min(sigma))

def get_Ainv(A: np.array):
    '''
    Calculate the inverse of matrix A, if it exists.

    Utilizing the SVD (U, Sigma, V) of matrix A, this function computes the
    inverse of A, raising an exception if A is singular (non-invertible).

    Parameters:
        A (np.array): The matrix to invert.

    Returns:
        A_inv (np.array): Inverse of matrix A.

    Raises:
        Exception: If A is singular (non-invertible).
    '''
    U, Sigma, V = svd(A)

    rows = A.shape[0]
    cols = A.shape[1]

    for i in range(min(cols, rows)):
        if Sigma[i][i] == 0:
            raise Exception(
                'Error: Input matrix is singular and cannot be inverted')
            

    if (np.diag(Sigma)).any() == 0:
        raise Exception('Error: Input matrix is singular and cannot be inverted') 
    else:
        sigma_inv = np.zeros((cols, rows))
        for i in range(min(cols, rows)):
                sigma_inv[i][i] = 1/Sigma[i][i]

    return V @ sigma_inv @ U.T
    
def pretty_matrix(A: np.array, name: str):
    '''
    Display a matrix or vector with a friendly format.

    This function takes a matrix or vector A and its name as inputs and 
    prints them to the console in a neat and readable format. The display 
    suppresses scientific notation and utilizes a precision of 5 decimal 
    places.

    Parameters:
        A (np.array): Matrix or vector to display.
        name (str): Descriptive name to print above the displayed matrix/vector.
    Returns:
        None
    '''
    np.set_printoptions(suppress=True, precision=5)
    if (np.ndim(A) == 1):
        print('-'*40, f'\n{name} Vector:\n')
    else:
        print('-'*40, f'\n{name} Matrix:\n')
    
    print(f'{A}\n')

def main():
    A = input_matrix()
    # A = [[1, 0], [1, 0]]  # Invertible Matrix
    # A = [[3, 2, 2], [2, 3, -2]]
    A = np.array(A)

    U, Sigma, V = svd(A)
    A_inv = get_Ainv(A)
    condition_num = get_condition(A)

    pretty_matrix(A, "Input")
    pretty_matrix(U, "U")
    pretty_matrix(Sigma, "Sigma")
    pretty_matrix(V.T, "V.T")
    pretty_matrix((U @ Sigma) @ V.T, "U @ Sigma @ V.T (should be A)")

    print('-'*40, f'\nCondition Number: ', ('%.5f' % condition_num))

    pretty_matrix(A_inv, "A Inverse")

    # Comparing to Blackbox
    U_b, Sigma_b, V_b = np.linalg.svd(A)
    pretty_matrix(A, "Blackbox Input")
    pretty_matrix(U_b, "Blackbox U")
    pretty_matrix(Sigma_b, "Blackbox Sigma")
    pretty_matrix(V_b, "Blackbox V")

# Example usage:
if __name__ == "__main__":
    main()
    
