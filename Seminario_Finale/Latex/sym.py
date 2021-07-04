import sympy as sy
a12, a13, a14, a21, a23, a24, a31, a32, a34, a41, a42, a43, x, y, z, w = sy.symbols('a12 a13 a14 a21 a23 a24 a31 a32 a34 a41 a42 a43 x y z w')
eq1 = sy.Eq(x + a12 * y + a13 * z + a14 * w, 1)
eq2 = sy.Eq(a21 * x + y + a23 * z + a24 * w, 1)
eq3 = sy.Eq(a31 * x + a32 * y + z + a34 * w, 1)
eq4 = sy.Eq(a41 * x + a42 * y + a43 * z + w, 1)
ans = sy.solve((eq1, eq2, eq3, eq4), (x, y, z, w))
print(sy.latex(ans))
