#!/usr/bin/env python3

import numpy as np
from math import sqrt

def input_matrix():
    '''
    Asks user for data, returns a np matrix
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
    Finds condition number of matrix using spectral norm
    '''
    eigVal_U, _ = np.linalg.eig(A@A.T)
    sigma = np.sqrt(eigVal_U)
    
    return (max(sigma)/min(sigma))

def get_Ainv(A):
    U, Sigma, V = svd(A)

    rows = A.shape[0]
    cols = A.shape[1]

    if (np.diag(Sigma)).any() == 0:
        raise Exception('Error: Input matrix is singular and cannot be inverted') 
    else:
        sigma_inv = np.zeros((cols, rows))
        for i in range(min(cols, rows)):
            sigma_inv[i][i] = 1/Sigma[i][i]

        return V @ sigma_inv @ U.T
    
def pretty_matrix(A: np.array, name: str):
    np.set_printoptions(suppress=True, precision=5)
    if (np.ndim(A) == 1):
        print('-'*40, f'\n{name} Vector:\n')
    else:
        print('-'*40, f'\n{name} Matrix:\n')
    
    print(f'{A}\n')

def main():
    # A = input_matrix()
    A = [[3, 2, 2], [2, 3, -2]]
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
    
