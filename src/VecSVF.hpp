#pragma once


struct Dioder
{
    CSmoothFilter vt,eta,is;
    double Vt,Eta,Is;
    Dioder()
    {
        vt.setCutoff(1/0.2);
        eta.setCutoff(1/0.2);
        is.setCutoff(1/0.2);
        Vt = 0.0253;
        Eta = 1.68;
        Is = .105;
    }
    Vec8f Diode(Vec8f x, double Vt = 0.0253,double eta = 1.68,double Is = .105)
    {    
        return Is * (exp(0.1*x/(eta*Vt))-1);
    }
    Vec8f Tick(Vec8f In)
    {
        double v = vt.process(Vt);
        if( fabs(v - Vt) < 1e-3) Vt  += noise.randint(-10,10)/1000.0;
        if(Vt < 0.001) Vt = 0.001;
        if(Vt > 0.2) Vt = 0.2;
        double e = eta.process(Eta);
        if(fabs(e - Eta) < 1e-3) Eta += noise.randint(-10,10)/1000.0;
        if(Eta < 1.0) Eta = 1.0;
        if(Eta > 2.0) Eta = 2.0;
        double i = is.process(Is);
        if(fabs(i - Is) < 1e-3) Is += noise.randint(-10,10)/10000.0;
        if(Is < .1) Is = .1;
        if(Is > 0.5) Is = 0.5;
        return Diode(In,Vt,Eta,Is);
    }
};

template<typename SIMD>
struct VecSVF
{    
    double fc,fs,q,K;
    SIMD lp,hp,bp,ubp,shelf,notch,apf,peak;
    SIMD z1,z2;
    Dioder diode;

    VecSVF(double Fs, double Fc, double Q) {
        fc = Fc;
        fs = Fs;
        q = Q;
        if(q == 0) q = 0.001;
        z1 = z2 = 0.0;
        K = 1;        
    }
    void setCutoff(float f) {
        if(f < 0) f = 0;                
        fc = 0.995*f;        
    }
    void setQ(float Q) {   
        if(q <= 0.0) q = 0.5;
        q = Q;
        
    }

    SIMD Tick(SIMD I, double A = 1, double X = 1, double Y = 1, double Z = 1, double OA=1, double minC = -1, double maxC = 1)
    {                
        float wd = 2*M_PI*fc;
        float T  = 1/fs;
        float temp = X*wd*T/2;
        float wa;

        if(fabs(temp) != M_PI/2)
            wa = (2/T)*tan(temp);
        else
            wa = (2/T)*tan(temp*0.995);
        
        float g  = wa*T/2;
        SIMD xn = A*I;
        float R  = fabs(Y)*1.0/(2*(q));
        if(R == 0) R = 0.001;
        
        //if(xn < minC) xn = minC;
        //if(xn > maxC) xn = maxC;
        xn = diode.Tick(xn);
        hp = (xn - (2*R+Z*g)*z1 - z2) / (1.0 + Z*2*R*g + Z*Z*g*g);        
        bp = g*hp + z1;        
        lp = g*bp + z2;
        
        // not sure these work right yet
        ubp = 2 * R * bp;
        // dont know exactly what K is it's not explained
        shelf = xn + 2*K*R*bp;
        notch = xn - 2*R*bp;
        apf   = xn - 4*R*bp;
        peak  = lp - hp;

        // delay feedback
        z1 = (g*hp + bp);
        z2 = (g*bp + lp);
                
        return lp;
    }    
};

using VecSVF4 = VecSVF<Vec4f>;
using VecSVF8 = VecSVF<Vec8f>;