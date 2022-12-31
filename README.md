# std-stk

# Synthesis Toolkit (Stk)
* LuaJIT
* Python
* GNU Octave

# TStk<T>
* float
* double
* complex
* autodiff

# stk::octave
* stk in octave for vectorizing experiments
  
```octave
wave = zeros(1,256);
v    = (0:255)/256;
for j=1:1000
  t = -1^j*(sin(2*pi*j*v)/j);
  wave += t;
end
wave = 0.5 - wave;
plot(wave);
pause();
```

  
# Vectorize
* Object can be vectorized as function of vector
* Object can be vectorized as objects per channel
* func(vector)
* func(matrix)
* for_each(samples,Func)
* Vector(Filter) = Filters[x]
* Matrix(Filter) = Filters(i,j)
* Tick(samples,FilterVector)
* Tick(samples,FilterMatrix)
