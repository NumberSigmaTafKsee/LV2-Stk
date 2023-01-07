require('se')
a = se.FloatMatrix(3,3)
t=os.clock()
for i=1,1000000 do
r = a*a
end
print(os.clock()-t)
