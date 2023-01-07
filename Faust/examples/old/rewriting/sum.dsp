
// How to emulate Faust's seq, par, sum.
// x(k) is assumed to yield the kth signal.

xseq(1,x)	= x(0);
xseq(n,x)	= xseq(n-1,x) : x(n-1);

xpar(1,x)	= x(0);
xpar(n,x)	= xpar(n-1,x) , x(n-1);

xsum(1,x)	= x(0);
xsum(n,x)	= xsum(n-1,x) + x(n-1);

// These are all very similar. Abstracting
// on the binary "accumulator" function, we
// get the familiar fold(-left) function:

fold(1,f,x)	= x(0);
fold(n,f,x)	= f(fold(n-1,f,x),x(n-1));

// Now seq, par, sum can be defined as:

fseq(n)		= fold(n,\(x,y).(x:y));
fpar(n)		= fold(n,\(x,y).(x,y));
fsum(n)		= fold(n,+);

// Often it is more convenient to specify
// parameters as a Faust tuple. We can match
// against the (x,y) pattern to decompose
// these.

vfold(f,(x,y))	= f(vfold(f,x),y);
vfold(f,x)	= x;

// Tuple version of seq, par, sum:

vseq		= vfold(\(x,y).(x:y));
vpar		= vfold(\(x,y).(x,y));
vsum		= vfold(+);

// Example: sum of sinusoids.

import("music.lib");

f0		= 440;
a(0)		= 1;
a(1)		= 0.5;
a(2)		= 0.3;
h(i)		= a(i)*osc((i+1)*f0);

v		= hslider("vol", 0.3, 0, 1, 0.01);

process		= v*fsum(3,h);
//process	= v*xsum(3,h);
//process	= v*vsum((h(0),h(1),h(2)));
