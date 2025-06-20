import numpy.polynomial.polynomial as poly
import scipy.signal

def round_to(num,digits): #digits = decimal digits
    return (round(num*(10**digits))/10**digits)
    
def process(str): #in form 5*x^4+1*x^0
    text = ''
    for i in str:
        if i != ' ':
            text+=i
    processed_ls = text.split('+')
    
    for i in range(len(processed_ls)):
        ls = [processed_ls[i].split('*')[0],'x',processed_ls[i].split('^')[1]]
        processed_ls[i] = ls
    
    return processed_ls

def anti_process(ls):  #in processed form
    normal_form = ""
    exp = ls.copy()
     #high_exp is the highest power of term in expression not included in normal form
    while len(exp)>0: #len(exp) has the number of terms
        #pick highest and include
        
        high_exp = max_power(exp)
        for j in range(len(exp)):
            if int(exp[j][2]) == high_exp:
                if normal_form != "":
                    normal_form += "+"
                normal_form += str(exp[j][0])+"*"+str(exp[j][1])+"^"+str(exp[j][2])
            exp.pop(j)
            break
    return normal_form

def display_mat(mat): #makes a matrix consisting of expressions into displayable format
    displayed_mat = []
    for i in range(len(mat)):
        row = []
        for j in range(len(mat[0])):
            row.append(anti_process(mat[i][j]))
        displayed_mat.append(row)
    return displayed_mat
            

def max_power(exp): #exp1,2 are in processed form ie [['1', 'x', '3'], ['3', 'x', '2'], ['3', 'x', '1'], ['1', 'x', '0']]
    highest = 0
    for i in range(len(exp)):
        if int(exp[i][2]) > highest:
            highest = int(exp[i][2])
    return highest

def create_zeropoly(exp1,exp2):
    l1 = max_power(exp1)
    l2 = max_power(exp2)
            
    max_exp = l1+l2
    poly = []
    for i in range(max_exp+1,-1,-1):
        poly.append(['0','x',str(i)])
    return poly

def fill(exp,max):  #exp1,2 are in processed form ie [['1', 'x', '3'], ['3', 'x', '2'], ['3', 'x', '1'], ['1', 'x', '0']] max is max power
    ls=exp
    powers = []
    for i in range(len(exp)):
        powers.append(int(exp[i][2]))
    temp = list(range(0,max+1))
    for i in temp:
        if i not in powers:
            ls.append(['0','x',str(i)])
    return ls

def transpose(mat):
    m = len(mat)
    n = len(mat[0])
    t = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(mat[j][i])
        t.append(row)
    return t
            
             
def polyadd(exp1,exp2): #exp1,2 are in processed form ie [['1', 'x', '3'], ['3', 'x', '2'], ['3', 'x', '1'], ['1', 'x', '0']]
    sum = []
    maximum = max([max_power(exp1),max_power(exp2)])
    e1 = fill(exp1,maximum)
    e2 = fill(exp2,maximum)
    for i in range(len(e1)):
        for j in range(len(e2)):
            if int(e1[i][2]) == int(e2[j][2]):
                sum.append([str(int(e1[i][0])+int(e2[j][0])),'x',e1[i][2]])               
    return sum
            
                        
def polyproduct(exp1,exp2): #mutiply 2 polynomials - exp1,2 are in processed form
    u1 = exp1
    u2 = exp2
    ans = create_zeropoly(u1,u2)
    for i in range(len(u1)):
        for j in range(len(u2)):
            product = [int(u1[i][0])*int(u2[j][0]),'x',int(u1[i][2])+int(u2[j][2])]
            for k in range(len(ans)):
                if int(ans[k][2]) == int(u1[i][2])+int(u2[j][2]):
                    ans[k][0] = str(int(u1[i][0])*int(u2[j][0])+int(ans[k][0]))
    return ans
            
            
def det_order2(mat):
    return polyadd((polyproduct(mat[0][0],mat[1][1])),polyproduct(polyproduct(mat[0][1],mat[1][0]),process("-1*x^0")))

def make_submat(mat,delrow,delcol):
    sub = []
    for i in range(len(mat)):
        if i != delrow:
            row = []
            for j in range(len(mat)):
                if j != delcol:
                    row.append(mat[i][j])
            sub.append(row)
    return sub

def det(mat): # mat is in processed form
    sub = []
    deter = [['0','x','0']]
    if len(mat)>2:
        for i in range(len(mat)):
            sub.append(make_submat(mat,0,i))
            deter = polyadd(polyproduct(polyproduct(mat[0][i],det(sub[i])),[[str(((((i+1)%2)*2)-1)),'x','0']]),deter)
        return deter
    elif len(mat) == 2:
        return det_order2(mat)
    else:
        print("newb")
        
def det_to_num(exp): #exp is output to the det function which is in processed form, used to find numerical det
    for term in exp:
        if term[2] == '0':
            return float(term[0])
     
def process_mat(inp):
    mat = []
    for i in range(len(inp)):
        row = []
        for j in range(len(inp[0])):
            row.append([[str(inp[i][j]),'x','0']])
        mat.append(row)
    return mat

def anti_process_mat(mat):
    num_mat = []
    for i in mat:
        row = []
        for j in i:
            row.append(det_to_num(j))
        num_mat.append(row)
    return num_mat
            
def identity(order):
    mat = []
    for i in range(order):
        row = []
        for j in range(order):
            if i == j:
                row.append([['-1','x','1']])
            else:
                row.append([['0','x','0']])
        mat.append(row)
    return mat

def mat_add(mat1,mat2):
    mat = []
    for i in range(len(mat1)):
        row = []
        for j in range(len(mat1[0])):
            row.append(polyadd(mat1[i][j],mat2[i][j]))
        mat.append(row)
    return mat
          
def characteristic_poly(mat):
    A = process_mat(mat)
    xI = identity(len(mat))
    poly = det(mat_add(A,xI))
    return poly

def coefficient_ls(exp):
    coefficient_mat = []
    high = max_power(exp)
    for i in range(high,-1,-1):
        for j in range(high+1):
            if int(exp[j][2]) == i and not(coefficient_mat==[] and int(exp[j][0]) == 0):
                coefficient_mat.append(int(exp[j][0]))             
    if coefficient_mat == []:
        coefficient_mat.append(0)    
    return coefficient_mat

def inv_coefficient(ls): #ls is the coeffiecient list in descending order of highest power
    exp = []
    for i in range(len(ls)):
        exp.append([str(ls[i]),'x',str(len(ls)-i-1)])
    return exp
        
def find_roots(coefficient_ls): #characterstic_ls [4,4,1] is 1x^2+4x+4
    roots = []
    root = poly.polyroots(coefficient_ls)
    for i in range(len(root)):
        if isinstance(root[i],complex): #isinstance function
            roots.append(complex(round_to(float(root[i].real),5),root[i].imag))
        else:
            roots.append((round_to(float(root[i].real),5)))
    roots = list(set(roots))
    return roots

def eigen_values(inp):
    values = find_roots(coefficient_ls(characteristic_poly(inp))[::-1]) 
    print("The Eigen Values are: ",end='')
    for i in values:
        print(i,end=' ')  
    print()
    return

def adjoint(mat): #input must be a processed matrix
    adj = []
    for i in range(len(mat)):
        row = []
        for j in range(len(mat[0])):
            row.append(polyproduct(det(make_submat(mat,i,j)),process(str((((i+j)%2)*-2)+1)+"*x^0")))
        adj.append(row)
    adj = transpose(adj)
    return adj

def inv(mat): #must be a processed matrix
    adj = adjoint(mat)
    d = det(mat)
    inv = []
    for i in range(len(mat)):
        row = []
        for j in range(len(mat[0])):
            divide_ls = list(scipy.signal.deconvolve(coefficient_ls(adj[i][j]),coefficient_ls(d))[0]) #consists of coefficients
            row.append(inv_coefficient(divide_ls))
        inv.append(row)
    return inv

def approximate(num,n): #3-n decimal place
    return int(num*10**n)/10**n

def approx_mat(mat,n): #integer entries
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            mat[i][j] = approximate(mat[i][j],n)
    return mat

def inp_mat(m,n): #where m and n is number of rows and columns 
    mat = []
    for i in range(m):
        row = []
        for j in range(n):
            row.append(process(input(f"item at row {i}, column {j} (in form 5*x^4+1*x^0) :  ")))
        mat.append(row)
    return mat #mat is in processed form

def print_mat(mat):
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            print(mat[i][j],end=" ")  
        print()

inp = [[1,2,3,6],[2,3,4,7],[4,5,9,8],[10,11,12,13]]#SQUARE MATRIX 
inp2 = [[0,-1],[1,0]] 
inp3 = [[2,3,4],[5,8,6],[0,4,2],[1,7,6]]
#eigen_values([[0,1],[-1,0]])
#eigen_values(inp2)
#print(characteristic_poly(inp2))
#print(process_mat(transpose(inp3)))
#print(det(process_mat(inp2)),"       ",inv_coefficient(coefficient_ls(det(process_mat(inp2)))))
#print(anti_process_mat(adjoint(process_mat(inp))))
#print(scipy.signal.fftconvolve([1,2,4,5],[2,7,6,2]))
#print(list(scipy.signal.deconvolve([2,11,28,52,63,38,10],[1,2,4,5])[0]))
#print(inv_coefficient([7,2,8,4]))
#print(anti_process_mat(adjoint(process_mat(inp))))
#print(det_to_num(det(process_mat(inp))))

"""
Features:
    Transpose
    det
    add
    characeristic poly
    eigen values
    adjoint 
    inverse
"""

features_available = ["Transpose","det","characteristic poly","eigen values","adjoint","inverse","add"]
number_of_matrix = int(input("unary operation or binary? Enter 1 or 2 :  "))
expressions = int(input("Matrix contains expressions? (1 for yes, 0 for no) :  "))

if number_of_matrix == 2:
    del features_available[0:6]
    row1 = int(input("number of rows of matrix 1:  "))
    column1 = int(input("number of columns of matrix 1:  "))
    row2 = int(input("number of rows of matrix 2:  "))
    column2 = int(input("number of columns of matrix 2:  "))
    print()
else:
    del features_available[6:]
    rows = int(input("number of rows of matrix:  "))
    columns = int(input("number of columns of matrix:  "))
    print()
    
    mat = inp_mat(rows,columns)
    if rows != columns:
        del features_available[1:]
        print("Transpose: \n")
        if expressions == 0:
            print_mat(transpose(mat))
        else:
            print_mat(display_mat(transpose(mat)))
        
    else:       
        if det(mat) == 0:
            features_available.remove("inverse")






