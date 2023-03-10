/*
 * Presence and Shelve filters as given in
 *   James A. Moorer
 *   The manifold joys of conformal mapping:
 *   applications to digital filtering in the studio
 *   JAES, Vol. 31, No. 11, 1983 November
 */

#define SPN MINDOUBLE

double bw2angle(a,bw)
double a,bw;
{
  double T,d,sn,cs,mag,delta,theta,tmp,a2,a4,asnd;

  T = tan(2.0*PI*bw);
  a2 = a*a;
  a4 = a2*a2;
  d = 2.0*a2*T;
  sn = (1.0 + a4)*T;
  cs = (1.0 - a4);
  mag = sqrt(sn*sn + cs*cs);
  d /= mag;
  delta = atan2(sn,cs);
  asnd = asin(d);
  theta = 0.5*(PI - asnd - delta);
  tmp = 0.5*(asnd-delta);
  if ((tmp > 0.0) && (tmp < theta)) theta = tmp;
  return(theta/(2.0*PI));
}

void presence(cf,boost,bw,a0,a1,a2,b1,b2)
double cf,boost,bw,*a0,*a1,*a2,*b1,*b2;
{
  double a,A,F,xfmbw,C,tmp,alphan,alphad,b0,recipb0,asq,F2,a2plus1,ma2plus1;

  a = tan(PI*(cf-0.25));
  asq = a*a;
  A = pow(10.0,boost/20.0);
  if ((boost < 6.0) && (boost > -6.0)) F = sqrt(A);
  else if (A > 1.0) F = A/sqrt(2.0);
  else F = A*sqrt(2.0);
  xfmbw = bw2angle(a,bw);

  C = 1.0/tan(2.0*PI*xfmbw);
  F2 = F*F;
  tmp = A*A - F2;
  if (fabs(tmp) <= SPN) alphad = C;
  else alphad = sqrt(C*C*(F2-1.0)/tmp);
  alphan = A*alphad;

  a2plus1 = 1.0 + asq;
  ma2plus1 = 1.0 - asq;
  *a0 = a2plus1 + alphan*ma2plus1;
  *a1 = 4.0*a;
  *a2 = a2plus1 - alphan*ma2plus1;

  b0 = a2plus1 + alphad*ma2plus1;
  *b2 = a2plus1 - alphad*ma2plus1;

  recipb0 = 1.0/b0;
  *a0 *= recipb0;
  *a1 *= recipb0;
  *a2 *= recipb0;
  *b1 = *a1;
  *b2 *= recipb0;
}

void shelve(cf,boost,a0,a1,a2,b1,b2)
double cf,boost,*a0,*a1,*a2,*b1,*b2;
{
  double a,A,F,tmp,b0,recipb0,asq,F2,gamma2,siggam2,gam2p1;
  double gamman,gammad,ta0,ta1,ta2,tb0,tb1,tb2,aa1,ab1;

  a = tan(PI*(cf-0.25));
  asq = a*a;
  A = pow(10.0,boost/20.0);
  if ((boost < 6.0) && (boost > -6.0)) F = sqrt(A);
  else if (A > 1.0) F = A/sqrt(2.0);
  else F = A*sqrt(2.0);

  F2 = F*F;
  tmp = A*A - F2;
  if (fabs(tmp) <= SPN) gammad = 1.0;
  else gammad = pow((F2-1.0)/tmp,0.25);
  gamman = sqrt(A)*gammad;

  gamma2 = gamman*gamman;
  gam2p1 = 1.0 + gamma2;
  siggam2 = 2.0*sqrt(2.0)/2.0*gamman;
  ta0 = gam2p1 + siggam2;
  ta1 = -2.0*(1.0 - gamma2);
  ta2 = gam2p1 - siggam2;

  gamma2 = gammad*gammad;
  gam2p1 = 1.0 + gamma2;
  siggam2 = 2.0*sqrt(2.0)/2.0*gammad;
  tb0 = gam2p1 + siggam2;
  tb1 = -2.0*(1.0 - gamma2);
  tb2 = gam2p1 - siggam2;

  aa1 = a*ta1;
  *a0 = ta0 + aa1 + asq*ta2;
  *a1 = 2.0*a*(ta0+ta2)+(1.0+asq)*ta1;
  *a2 = asq*ta0 + aa1 + ta2;

  ab1 = a*tb1;
  b0 = tb0 + ab1 + asq*tb2;
  *b1 = 2.0*a*(tb0+tb2)+(1.0+asq)*tb1;
  *b2 = asq*tb0 + ab1 + tb2;

  recipb0 = 1.0/b0;
  *a0 *= recipb0;
  *a1 *= recipb0;
  *a2 *= recipb0;
  *b1 *= recipb0;
  *b2 *= recipb0;
}

void initfilter(f)
filter *f;
{
  f->x1 = 0.0;
  f->x2 = 0.0;
  f->y1 = 0.0;
  f->y2 = 0.0;
  f->y = 0.0;
}

void setfilter_presence(f,freq,boost,bw)
filter *f;
double freq,boost,bw;
{
  presence(freq/(double)SR,boost,bw/(double)SR,
           &f->cx,&f->cx1,&f->cx2,&f->cy1,&f->cy2);
  f->cy1 = -f->cy1;
  f->cy2 = -f->cy2;
}

void setfilter_shelve(f,freq,boost)
filter *f;
double freq,boost;
{
  shelve(freq/(double)SR,boost,
   &f->cx,&f->cx1,&f->cx2,&f->cy1,&f->cy2);
  f->cy1 = -f->cy1;
  f->cy2 = -f->cy2;
}

void setfilter_shelvelowpass(f,freq,boost)
filter *f;
double freq,boost;
{
  double gain;

  gain = pow(10.0,boost/20.0);
  shelve(freq/(double)SR,boost,
   &f->cx,&f->cx1,&f->cx2,&f->cy1,&f->cy2);
  f->cx /= gain; 
  f->cx1 /= gain; 
  f->cx2 /= gain; 
  f->cy1 = -f->cy1;
  f->cy2 = -f->cy2;
}

/*
 * As in ''An introduction to digital filter theory'' by Julius O. Smith
 * and in Moore's book; I use the normalized version in Moore's book.
 */
void setfilter_2polebp(f,freq,R)
filter *f;
double freq,R;
{
  double theta;

  theta = 2.0*PI*freq/(double)SR;
  f->cx = 1.0-R;
  f->cx1 = 0.0;
  f->cx2 = -(1.0-R)*R;
  f->cy1 = 2.0*R*cos(theta);
  f->cy2 = -R*R;
}

/*
 * As in
 *   Stanley A. White
 *   Design of a digital biquadratic peaking or notch filter
 *   for digital audio equalization
 *   JAES, Vol. 34, No. 6, 1986 June
 */
void setfilter_peaknotch(f,freq,M,bw)
filter *f;
double freq,M,bw;
{
  double w0,p,om,ta,d;

  w0 = 2.0*PI*freq;
  if ((1.0/sqrt(2.0) < M) && (M < sqrt(2.0))) {
    fprintf(stderr,"peaknotch filter: 1/sqrt(2) < M < sqrt(2)\n");
    exit(-1);
  }
  if (M <= 1.0/sqrt(2.0)) p = sqrt(1.0-2.0*M*M);
  if (sqrt(2.0) <= M) p = sqrt(M*M-2.0);
  om = 2.0*PI*bw;
  ta = tan(om/((double)SR*2.0));
  d = p+ta;
  f->cx = (p+M*ta)/d;
  f->cx1 = -2.0*p*cos(w0/(double)SR)/d;
  f->cx2 = (p-M*ta)/d;
  f->cy1 = 2.0*p*cos(w0/(double)SR)/d;
  f->cy2 = -(p-ta)/d;
}

/*
 * Some JAES's article on ladder filter.
 * freq (Hz), gdb (dB), bw (Hz)
 */
void setfilter_peaknotch2(f,freq,gdb,bw)
filter *f;
double freq,gdb,bw;
{
  double k,w,bwr,abw,gain;

  k = pow(10.0,gdb/20.0);
  w = 2.0*PI*freq/(double)SR;
  bwr = 2.0*PI*bw/(double)SR;
  abw = (1.0-tan(bwr/2.0))/(1.0+tan(bwr/2.0));
  gain = 0.5*(1.0+k+abw-k*abw);
  f->cx = 1.0*gain;
  f->cx1 = gain*(-2.0*cos(w)*(1.0+abw))/(1.0+k+abw-k*abw);
  f->cx2 = gain*(abw+k*abw+1.0-k)/(abw-k*abw+1.0+k);
  f->cy1 = 2.0*cos(w)/(1.0+tan(bwr/2.0));
  f->cy2 = -abw;
}

double applyfilter(f,x)
filter *f;
double x;
{
  f->x = x;
  f->y = f->cx * f->x + f->cx1 * f->x1 + f->cx2 * f->x2
    + f->cy1 * f->y1 + f->cy2 * f->y2;
  f->x2 = f->x1;
  f->x1 = f->x;
  f->y2 = f->y1;
  f->y1 = f->y;
  return(f->y);
}
