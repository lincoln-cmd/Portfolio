# tutorial from: https://www.youtube.com/watch?v=PZJ5mjQyxR8

import numpy as np

'''
grid = [[5,3,0,0,7,0,0,0,0],
       [6,0,0,1,9,5,0,0,0],
       [0,9,8,0,0,0,0,6,0],
       [8,0,0,0,6,0,0,0,3],
       [4,0,0,8,0,3,0,0,1],
       [7,0,0,0,2,0,0,0,6],
       [0,6,0,0,0,0,2,8,0],
       [0,0,0,4,1,9,0,0,5],
       [0,0,0,0,8,0,0,0,0]]

'''

grid = []

for i in range(9):
    #n = list(map(int, input().split()))
    n = list()
    n2 = str(input())
    for j in n2:
        n.append(int(j))
    grid.append(n)


def possible(row, column, number):
    global grid
    # Is the number appearing in the given row?
    for i in range(0, 9):
        if grid[row][i] == number:
            return False
        
    
    # Is the number appearing in the given column?
    for i in range(0, 9):
        if grid[i][column] == number:
            return False
    
    # Is the number appearing in the given square?
    x0 = (column // 3) * 3
    y0 = (row // 3) * 3
    for i in range(0, 3):
        for j in range(0, 3):
            if grid[y0 + i][x0 + j] == number:
                return False
            
    return True

def solve():
    global grid
    for row in range(0, 9):
        for column in range(0, 9):
            if grid[row][column] == 0:
                for number in range(1, 10):
                    if possible(row, column, number):
                        grid[row][column] = number
                        solve()
                        grid[row][column] = 0
                        
                return
            
    print(np.matrix(grid))
    print('More possible solutions')
    
solve()