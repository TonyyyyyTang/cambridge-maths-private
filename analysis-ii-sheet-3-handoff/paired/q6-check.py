import sympy as s
x,y,t=s.symbols('x y t', real=True)
f=x**3/(x*x+y*y)
assert s.cancel(f.subs({x:t*x,y:t*y},simultaneous=True)-t*f)==0
val=lambda a,b: f.subs({x:a,y:b})
assert val(1,1)!=val(1,0)+val(0,1)
for curve,expected in [((t,t+t*t),s.Rational(1,2)),((t*t,t**3),0)]:
    assert s.limit(f.subs({x:curve[0],y:curve[1]},simultaneous=True)/t,t,0)==expected
print("Signed homogeneity, failure of additivity and derivatives along curved paths verified.")
