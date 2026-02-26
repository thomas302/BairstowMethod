# Modifications to the baistow method implementation found at
# https://github.com/PowerUpMasters/BairstowMethod/blob/master/bairstow.py
# -*- coding:utf-8 -*-

# Bairstow Method implemented by PowerUp Masters
# WhatsApp: +57 (311)-7103897 | +57 (311)6301489

import cmath
import random

Num = float | int | complex
Coeff = list[Num]

'''
   Bairstow's Method where:
      r = Initial guess
      s = Initial guess
     roots = Empty Array
   a = Coefficient's vector
   g = Polinomial's degree   
'''
def roots_quadratic(a:Coeff) -> tuple[Num, Num]:
    # Calculates roots using the quadratic formula for 
    # polynomial of degree two and this form:
    # f(x) = a_0 + a_1 * x + a_2 *x^2

    D:Num = (a[1]**2.0)-(4.0)*(a[2])*(a[0])
    X1:Num = (-a[1] - cmath.sqrt(D))/(2.0*a[2])
    X2:Num = (-a[1] + cmath.sqrt(D))/(2.0*a[2])

    return X1, X2

def bairstow(a:Coeff, r:Num, s:Num, g:int, roots:list[Num]) -> None:
    if(g<1):
        return None

    if((g==1) and (a[1] != 0)):
        # Single root for equation of form:
        # f(x) = a_0 + a_1 * x
        roots.append((-a[0])/(a[1]))
        return None

    if(g==2):
        # Calculates roots using the quadratic formula for 
        # polynomial of degree two and this form:
        # f(x) = a_0 + a_1 * x + a_2 *x^2
        _x = roots_quadratic(a)
        roots.append(_x[0])
        roots.append(_x[1])
        return None

    n:int = len(a)
    b:Coeff = [0.0]*len(a)
    c:Coeff = [0.0]*len(a)

    b[n-1] = a[n-1]
    b[n-2] = a[n-2] + r*b[n-1]
    
    i = n - 3
    while(i>=0):
        b[i] = a[i] + r*b[i+1] + s*b[i+2]
        i = i - 1
    c[n-1] = b[n-1]
    c[n-2] = b[n-2] + r*c[n-1]
    i = n - 3
    while(i>=0):
        c[i] = b[i] + r*c[i+1] + s*c[i+2]
        i = i - 1

    Din:Num = ((c[2]*c[2])-(c[3]*c[1]))**(-1.0)

    r:Num = r + (Din)*((c[2])*(-b[1])+(-c[3])*(-b[0]))
    s:Num = s + (Din)*((-c[1])*(-b[1])+(c[2])*(-b[0]))

    if(abs(b[0])>1E-14 or abs(b[1])>1E-14):
        return bairstow(a,r,s,g,roots)
    if (g>=3):
        Dis = ((-r)**(2.0))-((4.0)*(1.0)*(-s))
        X1 = (r - (cmath.sqrt(Dis)))/(2.0)
        X2 = (r + (cmath.sqrt(Dis)))/(2.0)
        roots.append(X1)
        roots.append(X2)
        return bairstow(b[2:],r,s,g-2,roots)    


if __name__ == "__main__":
    g = int(input("\nEnter the rank of your polynomial: \n"))

    a_in = input("\nEnter Coefficients of polynomial seperated by commas: \n")
    roots:Coeff = []

    a:Coeff = [complex(x) for x in a_in.strip().split(",")]

    # --- Optional r ---
    r_in = input("\nEnter initial guess r (press Enter for random): ").strip()
    r: Num = complex(r_in) if r_in else random.random()

    # --- Optional s ---
    s_in = input("Enter initial guess s (press Enter for random): ").strip()
    s: Num = complex(s_in) if s_in else random.random()
    bairstow(a,r,s,g,roots)

    print("\nFound Roots => \n")
    for k, r in enumerate(roots):
        print("R" + str(k) + " = " + str(r))
